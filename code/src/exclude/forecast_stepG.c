/*
 * forecastG.c
 *
 *  Created on: May 17, 2012
 *      Author: camcat
 */

#include "forecast_stepG.h"


#define Max 100000
#define MaxTS 1000
#define Maxeq 10000

int forecast_stepG2_new(struct catalog cat, double *times, double **cmpdata, struct pscmp *DCFS, double tt0, double tt1, double Asig, double ta,
			int points[], double *out_NeX, double *NeT, double *Rate_end, int N, int NTS, int Neqks, double *gamma_init, double *R, int last){
//assumes stress grows linearly during each time step.
//NB assumes that times, cat, DCFS, NTS are always the same (each time function is called).
// if points==NULL, use sequence [1,2,3,...,N].
// cmbdata can be NULL, and will be ignored.
	//fixme check: indices of R[0...cat.Z-1] or {[1...cat.Z]

  double  tau, dtau_dt, dtau_dt0;
  int     TS0, TS1, n;
  double  Tpre=0, dt1;
  double  gamma, gamma0;
  static double **events;	//will contain times, magnitudes of all earthquakes (both from catalog and from DCFS).
  int Neq, err=0;
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
  double  epsilongamma=1e-200;	//minimum tau for evolving gamma (otherwise use gamma=gamma_0+dt/Asig); minimum gamma for calculating time step using frac*gamma*Asig.
  double frac=0.1;
  double a,b;
  double dt_tiny0=ta*1e-14; //ta*1e-14;	//frac:maximum fractional change in gamma allowed (used to compute tiny time steps when necessary). dt_tiny0 is value used when gamma0=0 to avoid infinite loop (or very small steps if exp(-t_tiny/t_a)=1, due to representation of doubles). //todo experiment with this!
  double dt_tiny_cum, dt_tiny;
  int exit_loop=0, reach_end;
  int nthreads;
  static int **indices;	//containes indices of events referred to DCFS and cat(with offset of 1 - see below).

  //todo could use analytical solution for time step if afterslip==0.
  //todo make function works with times2=NULL (for case afterslip==0).
  //todo simplify calculations for times before stresses.

  if (firsttimein==1){
		firsttimein=0;
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
			if (k<0 && events[1][i]>=tt0 && verbose_level>0) printf("**Warning: no time steps available before earthquake no. %d ** (forecast_stepG2_new.c)\n",i);
		}
  }

  if (Asig==0 && ta==0.0) return(0);	//in this case, function has only been called to setup variables above.


  for (int z=1; z<=cat.Z; z++) if (cat.t[z]>=tt0 && cat.t[z]<tt1) R[z]=0.0;
  NeX= (out_NeX)? out_NeX : NeXdum;

  // find time steps before each contributing earthquake;
  //if ((old_DCFS_addr!=new_DCFS_addr) && verbose_level>1) printf("**Warning: adress of DCFShas changed in forecast_stepG2_new.c\n",i);
  old_DCFS_addr=new_DCFS_addr;
  if (Rate_end) for(int m=1;m<=N;m++) ReX[m]=0.0;
  if (NeX) for(int m=1;m<=N;m++) NeX[m]=0.0;

  if (tt1<tt0 && verbose_level>1) printf("\n*** Warning: tt1<tt0 in forecast_stepG2_new.c  ***\n");
  if (times[NTS]>tt0 && times[NTS]<tt1 && verbose_level>1) {
	  printf("\n** Warning: time steps in forecast_stepG don't cover entire forecast range!**\n");
  }

  if (times[0]>=tt1){
	  if (NeT!=(double *) 0) *NeT=N*(tt1-tt0);
	  if (NeX) for(int m=1;m<=N;m++) NeX[m]=(tt1-tt0);
	  return (1);
 	  }

//TS0= First time step after tt0.
  TS0=0;
  while(TS0<NTS-1 && times[TS0]<=tt0) TS0++;
  if(times[TS0]<=tt0 && verbose_level>1) {
	  printf("\n*** Warning: times[TS0]<=tt0 in forecast_stepG2_new.c  ***\n");
  }

//TS1= Last time step before tt1.
  TS1=1;
  while(TS1<NTS-1 && times[TS1]<tt1) TS1++;
  if(times[TS1]>=tt1) TS1--;
  dt1=tt1-times[TS1];

  if (NeT!= (double *) 0) *NeT=0;
  if (Rate_end!= (double *) 0) *Rate_end=0;

//  printf ("max no of threads=%d.\n", omp_get_max_threads());

  Rprivate=dmatrix(0,omp_get_max_threads()-1, 0, cat.Z);
  for (int t=0; t<omp_get_max_threads(); t++){
	  for (int eq=0; eq<=cat.Z; eq++) Rprivate[t][eq]=0.0;
  }

#pragma omp parallel for private(err, n, gamma,gamma0,tau,j0, next_eqk, next_TS, t_now, t_pre, dt_tiny_cum, exit_loop, dt_tiny, reach_end, t_endstep, cat_i, DCFS_i, a, b, dtau_dt, dtau_dt0, counter_eqk)
  for(int m=1;m<=N;m++){

	nthreads=omp_get_num_threads();
//	printf("nthreads=%d\n",nthreads);
//	printf("here's thread no. %d.\n",omp_get_thread_num());
	//printf("Entered parallel section - %d threads running. \n",omp_get_num_threads());

	n=(points==0)? m : points[m];
	if (NeX) NeX[m]=Tpre;
	gamma=gamma_init[m];

	gamma0=gamma;
    j0=TS0;
    counter_eqk=0;
    t_now=tt0;
    reach_end=0;

    err=0;
    dtau_dt=dtau_dt0=Asig/ta;

    //new method (uses indices of earthquakes for specific gridpoint):

    while (counter_eqk<num_eqk[n] && events[1][which_eqk[n][counter_eqk]]<tt0) counter_eqk+=1;

    while (tol0<(tt1-t_now)){

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

		dt_tiny_cum=0.0;
		exit_loop=0;
		while (exit_loop==0){
			gamma0=gamma;
			dt_tiny= (gamma0>epsilongamma)? frac*gamma0*Asig : dt_tiny0;		//for gamma0<<1, this time step size gives |gamma -gamma0|=frac*gamma0;
			if (dt_tiny==0) dt_tiny=dt_tiny0;	//line above should prevent infinite loop, but this adds one more check.
			if (dt_tiny_cum+dt_tiny>=t_pre) {
				dt_tiny= t_pre-dt_tiny_cum;
				exit_loop=1;
			}
			tau=(cmpdata && j0-1>=0)? (Asig/ta)*dt_tiny+(dt_tiny/dt[j0-1])*cmpdata[j0-1][n] : (Asig/ta)*dt_tiny;
			gamma=(fabs(tau/Asig)>1e-10)? (gamma-dt_tiny/(tau))*exp(-tau/Asig)+dt_tiny/(tau) : gamma*(1-tau/Asig)+dt_tiny/Asig;
			if (NeX) NeX[m]+=(ta/Asig)*dt_tiny*0.5*(1/gamma+1/gamma0);
			dt_tiny_cum+=dt_tiny;
		}

		t_now+=t_pre;

		for(int j=j0;j<next_TS;j++){

			dt_tiny_cum=0.0;
			exit_loop=(t_pre==0);
			while (exit_loop==0){
				gamma0=gamma;
				dt_tiny= (gamma0>epsilongamma)? frac*gamma0*Asig : dt_tiny0;		//for gamma0<<1, this time step size gives |gamma -gamma0|=frac*gamma0;
				if (dt_tiny==0) dt_tiny=dt_tiny0;	//line above should prevent infinite loop, but this adds one more check.
				if (dt_tiny_cum+dt_tiny>=dt[j]) {
					dt_tiny= dt[j]-dt_tiny_cum;
					exit_loop=1;
				}
				tau=(cmpdata)? (Asig/ta)*dt_tiny+(dt_tiny/dt[j])*cmpdata[j][n] : (Asig/ta)*dt_tiny;
				gamma=(fabs(tau/Asig)>1e-10)? (gamma-dt_tiny/(tau))*exp(-tau/Asig)+dt_tiny/(tau) : gamma*(1-tau/Asig)+dt_tiny/Asig;
				if (NeX) NeX[m]+=(ta/Asig)*dt_tiny*0.5*(1/gamma+1/gamma0);
				dt_tiny_cum+=dt_tiny;
			}
			t_now+=dt[j];
		}

		t_pre=t_endstep-t_now;

		dt_tiny_cum=0.0;
		exit_loop=(t_pre<tol0);
		while (exit_loop==0){
			gamma0=gamma;
			dt_tiny= (gamma0>epsilongamma)? frac*gamma0*Asig : dt_tiny0;		//for gamma0<<1, this time step size gives |gamma -gamma0|=frac*gamma0;
			if (dt_tiny==0) dt_tiny=dt_tiny0;	//line above should prevent infinite loop, but this adds one more check.
			if (dt_tiny_cum+dt_tiny>=t_pre) {
				dt_tiny= t_pre-dt_tiny_cum;
				exit_loop=1;
			}
			tau=(cmpdata)? (Asig/ta)*dt_tiny+(dt_tiny/dt[next_TS])*cmpdata[next_TS][n] : (Asig/ta)*dt_tiny;
			gamma=(fabs(tau/Asig)>1e-10)? (gamma-dt_tiny/(tau))*exp(-tau/Asig)+dt_tiny/(tau) : gamma*(1-tau/Asig)+dt_tiny/Asig;
			if (NeX) NeX[m]+=(ta/Asig)*dt_tiny*0.5*(1/gamma+1/gamma0);
			dt_tiny_cum+=dt_tiny;
		}
		t_now+=t_pre;

		if (reach_end==0){
			gamma0=gamma;
			if (cat_i!=0) {
				Rprivate[omp_get_thread_num()][cat_i]+= cat_weights[n][counter_eqk]*(ta/Asig)/gamma;
			}
			gamma= (DCFS_i==-1 | DCFS_whichpt[n][counter_eqk]==0)? gamma : gamma*exp(-DCFS[DCFS_i].cmb[DCFS_whichpt[n][counter_eqk]]/Asig);
			if (isinf(gamma)==1){
				if (verbose_level>0) printf("*Warning: gamma==Inf, must choose larger Asig!*\n");
				err=1;
			}
		}
		j0=next_TS+1;
		counter_eqk+=1;
	}

	if (last) gamma_init[m]=gamma;
	if (Rate_end) ReX[m]=(ta/Asig)/gamma;

  }

  for (int t=0; t<nthreads; t++){
	//for (int eq=0; eq<Neq; eq++) R[indices[1][eq]+1]+=Rprivate[t][eq];
	  for (int eq=0; eq<=cat.Z; eq++) R[eq]+=Rprivate[t][eq];
  }
  if (Rate_end) *Rate_end=0.0;
  if (NeT) *NeT=0.0;
  for(int m=1;m<=N;m++){
	  if (NeT) *NeT+=NeX[m];
	  if (Rate_end) *Rate_end+=ReX[m];
  }

  free_dmatrix(Rprivate,0,omp_get_max_threads()-1, 0, Neq-1);
  return(err);

}
