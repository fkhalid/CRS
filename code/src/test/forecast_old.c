/* This is the same as forecast_stepG.c, but uses steps insted of linear evolution b/w time steps.
 * forecastG.c
 *
 *  Created on: May 17, 2012
 *      Author: camcat
 */

#include "forecast_old.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

#define Max 100000
#define MaxTS 1000
#define Maxeq 10000

int forecast_stepG2_old(struct catalog cat, double *times, double **cmpdata, struct pscmp *DCFS, double tt0, double tt1, double Asig, double ta,
			int points[], double *out_NeX, double *NeT, double *Rate_end, int N, int NTS, int Neqks, double *gamma_init, double *back_rate, double *R, int last){
//assumes stress grows as a step in the center of each time step.
//NB assumes that times, cat, DCFS, NTS are always the same (each time function is called).
// if points==NULL, use sequence [1,2,3,...,N].
// cmbdata can be NULL, and will be ignored.
// if backrate==1, will assume it's 1 for all points. If not, it should be an array containing NgridT times the ratio between avg rate and rate at that point. (so that the sum of back_rate is NgridT).
	//fixme check: indices of R[0...cat.Z-1] or {[1...cat.Z]

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

  double  tau, dtau_dt, dtau_dt00=Asig/ta, ta1;
  int     TS0, TS1, n;
  double  Tpre=0, dt1;
  double  gamma, gamma0, back_rate_n;
  static double **events;	//will contain times, magnitudes of all earthquakes (both from catalog and from DCFS).
  int Neq, err=0, errtot=0;
  int is_incat, is_inDCFS;
  int j0, next_eqk, counter_eqk;
  int next_TS, cat_i, DCFS_i;
  double *NeX;
  double **Rprivate;
  static double *dt, *NeXdum, *ReX;
  static int *TS_eqk;
  static int firsttimein=1;
  static int **which_eqk;	//contains list of eqk for each grid point, and index of the point within DCFS;
  static int *num_eqk;
  static float ** cat_weights;
  static int ** DCFS_whichpt;
  double t_pre, t_now, t_endstep;
  static long old_DCFS_addr=0;
  long new_DCFS_addr= (long) DCFS;
  double a,b;
  int reach_end;
  int nthreads, nthreadstot=omp_get_max_threads();
  static int **indices;	//contains indices of events referred to DCFS and cat(with offset of 1 - see below).
  int warning_printed=0;

  if (firsttimein==1){
		firsttimein=0;
		//todo: use previous information to match DCFS and cat.
		events=union_cats(cat.t+1, timesfrompscmp(DCFS,Neqks), cat.mag+1, magsfrompscmp(DCFS,Neqks), cat.Z, Neqks, 0.001, 0.3, &indices, &Neq);	//todo parameters should not be hardwired.
		NeXdum=dvector(1,N);
		ReX=dvector(1,N);

		//list of all events, total no. of events affecting each point (Neq: from combining DCFS and cat).
		which_eqk=imatrix(1,N,0,Neq-1);
		num_eqk=ivector(1,N);
		cat_weights=matrix(1,N,0,Neq-1);
		DCFS_whichpt=imatrix(1,N,0,Neq-1);
		dt=dvector(0,NTS-1);
		TS_eqk=ivector(0,Neq-1);

		for (int n=1; n<=N; n++) num_eqk[n]=0;

		for(int eq=0;eq<Neq;eq++){
			int counter1, counter2;
			counter1=counter2=1;
			cat_i=indices[1][eq]+1;	//	NB:cat_i[i]==0 means no events selected.
			DCFS_i=indices[2][eq];  //	NB:DCFS_i[i]==-1 means no events selected.
			for (int n=1; n<=N; n++){
				if (counter2 > DCFS[DCFS_i].nsel && counter1 > cat.ngrid[cat_i]) break;	//all points have already been included
				else{
					if (cat_i!=0) while (counter1 <=cat.ngrid[cat_i] && cat.ngridpoints[cat_i][counter1]<n) counter1+=1;
					if (DCFS_i!=-1) while (counter2 <=DCFS[DCFS_i].nsel && DCFS[DCFS_i].which_pts[counter2]<n) counter2+=1;
					is_incat= (cat_i!=0 && counter1 <= cat.ngrid[cat_i] && cat.ngridpoints[cat_i][counter1]==n) ? 1 : 0;
					is_inDCFS= (DCFS_i!=-1 && counter2 <= DCFS[DCFS_i].nsel && DCFS[DCFS_i].which_pts[counter2]==n) ? 1 : 0;
					if (is_incat | is_inDCFS){
						which_eqk[n][num_eqk[n]]=eq;
						cat_weights[n][num_eqk[n]]= (is_incat)? cat.weights[cat_i][counter1] : 0.0;
						DCFS_whichpt[n][num_eqk[n]]= (is_inDCFS)? counter2 : 0;
						num_eqk[n]+=1;
					}
				}
			}
		}

		for (int g=0; g<NTS; g++) dt[g]=times[g+1]-times[g];	//should be <=NTS?
		for (int i=0; i<Neq; i++){
			int k=0;
			while(k<NTS && times[k]<events[1][i]) k++;
			if(times[k]>=events[1][i]) k--;
			TS_eqk[i]=k;
			if(procId == 0) {
				if (k<0 && events[1][i]>=tt0) {
					print_screen("**Warning: no time steps available before earthquake no. %d ** (forecast_stepG2_new.c)\n",i);
					print_logfile("**Warning: no time steps available before earthquake no. %d ** (forecast_stepG2_new.c)\n",i);
				}
			}
		}
  }

  if (Asig==0 && ta==0.0) return(0);	//in this case, function has only been called to setup variables above.


  //for (int z=1; z<=cat.Z; z++) if (cat.t[z]>=tt0 && cat.t[z]<tt1) R[z]=0.0;	//fixme check that vector is always initialized beforehand.
  NeX= (out_NeX)? out_NeX : NeXdum;

  // find time steps before each contributing earthquake;
  old_DCFS_addr=new_DCFS_addr;
  if (Rate_end) for(int m=1;m<=N;m++) ReX[m]=0.0;
  if (NeX) for(int m=1;m<=N;m++) NeX[m]=0.0;

  if(procId == 0) {
	  if (tt1<tt0) {
		  print_screen("\n*** Warning: tt1<tt0 in forecast_stepG2_new.c  ***\n");
		  print_logfile("\n*** Warning: tt1<tt0 in forecast_stepG2_new.c  ***\n");
	  }
	  if (times[NTS]>tt0 && times[NTS]<tt1) {
		  print_screen("\n** Warning: time steps in forecast_stepG don't cover entire forecast range!**\n");
		  print_logfile("\n** Warning: time steps in forecast_stepG don't cover entire forecast range!**\n");
	  }
  }

// todo [coverage] this block is never tested
  if (times[0]>=tt1){
	  if (NeT!=(double *) 0) *NeT=N*(tt1-tt0);
	  if (NeX) for(int m=1;m<=N;m++) NeX[m]=(tt1-tt0);

	  return (0);
  }

//TS0= First time step after tt0.
  TS0=0;
  while(TS0<NTS-1 && times[TS0]<=tt0) TS0++;
  if(procId == 0) {
	  if(times[TS0]<=tt0 & !warning_printed) {
		  warning_printed=1;
		  print_screen("\n*** Warning: times[TS0]<=tt0 in forecast_stepG2_new.c  ***\n");
		  print_logfile("\n*** Warning: times[TS0]<=tt0 in forecast_stepG2_new.c  ***\n");
	  }
  }

//TS1= Last time step before tt1.
  TS1=1;
  while(TS1<NTS-1 && times[TS1]<tt1) TS1++;
  if(times[TS1]>=tt1) TS1--;
  dt1=tt1-times[TS1];

  if (NeT) *NeT=0;
  if (Rate_end) *Rate_end=0;


  Rprivate=dmatrix(0,nthreadstot-1, 0, cat.Z);
  for (int t=0; t<omp_get_max_threads(); t++){
	  for (int eq=0; eq<=cat.Z; eq++) Rprivate[t][eq]=0.0;
  }

  err=0;
#pragma omp parallel for firstprivate(err) private(n, gamma,gamma0,tau,j0, next_eqk, next_TS, t_now, t_pre, reach_end, t_endstep, cat_i, DCFS_i, a, b, dtau_dt, counter_eqk, back_rate_n, ta1) reduction(+:errtot)
  for(int m=1;m<=N;m++){

	nthreads=omp_get_num_threads();

	if (err!=0) continue;
	n=(points==0)? m : points[m];

	back_rate_n= (back_rate) ? back_rate[n] : 1.0;
	if (NeX) NeX[m]=Tpre;
	gamma=gamma_init[m];

	gamma0=gamma;
    j0=TS0;
    counter_eqk=0;
    t_now=tt0;
    reach_end=0;

    //err=0;
    dtau_dt=Asig/ta;

    //new method (uses indices of earthquakes for specific gridpoint):

    while (counter_eqk<num_eqk[n] && events[1][which_eqk[n][counter_eqk]]<tt0) counter_eqk+=1;

    while (t_now<tt1){

    	if (err!=0) break;	//error in a previous loop;
    	if (counter_eqk>=num_eqk[n]) reach_end=1;
    	else {
        	next_eqk=which_eqk[n][counter_eqk];
    		if (events[1][next_eqk]>=tt1) reach_end=1;
    	}

    	next_TS= reach_end ? TS1 : TS_eqk[next_eqk];
    	t_endstep= reach_end ? tt1 : events[1][next_eqk];
    	cat_i= reach_end ? 0 : indices[1][next_eqk]+1;
    	DCFS_i= reach_end ? 0 : indices[2][next_eqk];

    	//find rates up to next barrier (which is the smallest between next earthquake, next time step or tt1).
    	t_pre= fmin(t_endstep, times[j0]) - t_now;

    	dtau_dt=(cmpdata && j0-1>=0 )? (Asig/ta)+cmpdata[j0-1][n]/dt[j0-1] : (Asig/ta);

    	if (t_pre>0.0){
			tau=dtau_dt*t_pre;
			gamma=(gamma+t_pre/(2.0*Asig))*exp(-tau/Asig)+t_pre/(2.0*Asig);
    	}

    	t_now+=t_pre;

		for(int j=j0;j<next_TS;j++){
			dtau_dt=(cmpdata)? (Asig/ta)+cmpdata[j][n]/dt[j] : (Asig/ta);
			tau=dtau_dt*dt[j];
			gamma=(gamma+dt[j]/(2.0*Asig))*exp(-tau/Asig)+dt[j]/(2.0*Asig);
			t_now+=dt[j];
		}

		t_pre=t_endstep-t_now;

		if (t_pre>0.0) {
			dtau_dt=(cmpdata)? (Asig/ta)+(1.0/dt[next_TS])*cmpdata[next_TS][n] : (Asig/ta);
			tau=dtau_dt*t_pre;
			gamma= (gamma+t_pre/(2.0*Asig))*exp(-tau/Asig)+t_pre/(2.0*Asig);
		}

		t_now+=t_pre;

		if (reach_end==0){
			gamma0=gamma;
			if (cat_i!=0) {
				Rprivate[omp_get_thread_num()][cat_i]+= back_rate_n*cat_weights[n][counter_eqk]*(ta/Asig)/gamma;
			}
			gamma= (DCFS_i==-1 | DCFS_whichpt[n][counter_eqk]==0)? gamma : gamma*exp(-DCFS[DCFS_i].cmb[DCFS_whichpt[n][counter_eqk]]/Asig);
			if (isinf(gamma)){
				print_screen("*Warning: gamma==Inf, must choose larger Asig!*\n");
				print_logfile("*Warning: gamma==Inf, must choose larger Asig!*\n");
				err=1;
				errtot+=1;
			}
		}
		j0=next_TS+1;
		counter_eqk+=1;
	}

	if (last) gamma_init[m]=gamma;
	if (Rate_end) ReX[m]=back_rate_n*(ta/Asig)/gamma;

  }

  if (R) for (int t=0; t<nthreads; t++){
	  for (int eq=0; eq<=cat.Z; eq++) R[eq]+=Rprivate[t][eq];
  }
  if (Rate_end) *Rate_end=0.0;
  if (NeT) *NeT=0.0;
  for(int m=1;m<=N;m++){
	  if (NeT) *NeT+=NeX[m];
	  if (Rate_end) *Rate_end+=ReX[m];
  }

  free_dmatrix(Rprivate,0,nthreadstot-1, 0, cat.Z);
  return(errtot);

}

