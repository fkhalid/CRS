/*
 * ETASforecast.c
 *
 *  Created on: May 17, 2012
 *      Author: camcat
 */

#include "ETASforecast.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>

#include "../defines.h"
#include "../inp_out/write_csep_forecast.h"
#include "../util/nrutil.h"
#include "ETASspatialkernel.h"
#include "struct_conversions.h"

#define Max 100000
#define MaxTS 1000
#define Maxeq 10000

int ETASforecast (struct pscmp *DCFS, struct eqkfm *eqkfm, struct crust crst, int NTSdisc, int NgridT, int Nm, struct catalog cat,
		double *tts, int Ntts, double tw, char *print_forex0, char *print_foret, char *print_LL, double etas_p,  double etas_c,  double etas_K,  double etas_d,
		double etas_q,  double etas_mu, double etas_alpha, double Mc){

	char fname[120];
	char print_forex[120];
	static double *dumrate, *rate, *ev_x=0;
	double sum, sum1, sum2, integral;
	double Ldum;
	double fin_rate;
	double *nev, *rev, *ev_x_dum;
	int N, NgridT_out= crst.uniform ? (crst.nLat_out*crst.nLon_out*crst.nD_out) : NgridT;
	int err=0;
	int current_main;
	double tnow, tt0, tt1;
	FILE *fLLev;
	FILE *fforet_avg;


	dumrate=dvector(1,cat.Z);
	rate=dvector(1,cat.Z);
	ev_x=dvector(1,NgridT);
	ev_x_dum=dvector(1,NgridT);
	nev=dvector(1,Ntts);
	rev=dvector(1,Ntts);

	for (int t=1; t<=Ntts; t++) {
		nev[t]=0.0;
		rev[t]=0.0;
	}
	N=0;

	tt0=tts[0];
	tt1=tts[Ntts];
	for(int i=1;i<=cat.Z;i++) if(cat.t[i]>=tt0 && cat.t[i]<=tt1) N+=1;

	sum=sum1=sum2=0.0;
	integral=0.0;
	for(int i=1;i<=cat.Z;i++) dumrate[i]=rate[i]=0.0;

	sum=sum1=sum2=0;
	err=0;
	if (print_forex0) {
		sprintf(print_forex,"%s.dat", print_forex0);
	}
	if (print_foret) {
		sprintf(fname, "%s.dat",print_foret);
		fforet_avg=fopen(fname,"w");
	}

	if (print_LL) {
		sprintf(fname, "%s.dat",print_LL);
		fLLev=fopen(fname,"w");
	}

	ETASspatialkernel(DCFS, eqkfm, crst, NTSdisc, tts[0], tts[Ntts], etas_d, etas_q, etas_alpha, Mc);

	for (int n=1; n<=NgridT; n++) ev_x[n]=0.0;
	for(int i=1;i<=cat.Z;i++) dumrate[i]=0.0;

	fflush(stdout);

	//todo delete
	forecastETAS(cat, DCFS, -245.0, 0.0, 0, ev_x_dum, &sum, &fin_rate, NgridT, NTSdisc, dumrate, 0.00972, 0.95, 0.052894);

	tt0=tts[0];
	tt1=tts[Ntts];

	current_main=0;
	tnow=tts[0];
		for (int t=1; t<=Ntts; t++){
		//Calculate seismicity evolution (skipping a time window after each mainshock):
		tt0=tts[t-1];
		tt1=tts[t];
		while (current_main<Nm && eqkfm[current_main].t<tt1){
			if (tnow<eqkfm[current_main].t){
				forecastETAS(cat, DCFS, tnow, eqkfm[current_main].t, 0, ev_x_dum, &sum, &fin_rate, NgridT, NTSdisc, dumrate, etas_c, etas_p, etas_K);
				for (int i=1; i<=NgridT; i++) ev_x[i]+=ev_x_dum[i];
				nev[t]=sum;
				rev[t]=fin_rate;
				forecastETAS(cat, DCFS, eqkfm[current_main].t, eqkfm[current_main].t+tw, 0, ev_x_dum, &sum, &fin_rate, NgridT, NTSdisc, dumrate, etas_c, etas_p,etas_K);
				tnow=eqkfm[current_main].t+tw;
			}
			else if (tnow<eqkfm[current_main].t+tw){
				forecastETAS(cat, DCFS, tnow, eqkfm[current_main].t+tw, 0, ev_x_dum, &sum, &fin_rate, NgridT, NTSdisc, dumrate, etas_c,etas_p, etas_K);
				tnow=eqkfm[current_main].t+tw;
			}
			current_main+=1;
			if (err) break;
		}
		if (tnow<tt1){
			forecastETAS(cat, DCFS, tnow, tt1, 0, ev_x_dum, &sum, &fin_rate, NgridT, NTSdisc, dumrate, etas_c, etas_p, etas_K);
			for (int i=1; i<=NgridT; i++) ev_x[i]+=ev_x_dum[i];
			nev[t]=sum;
			rev[t]=fin_rate;
			tnow=tt1;
		}
		for(int i=1;i<=cat.Z;i++) if(cat.t[i]>=tt0 && cat.t[i]<tt1) rate[i]+=dumrate[i];
	}

	//calculate average rate and LL:
	if (!err){

		if (print_foret) {
			for (int t=1; t<=Ntts; t++) {
				if (tts[t-1]<tts[0]+tw) fprintf(fforet_avg, "%lf\t%lf\t%lf\n",tts[t],0.0, 0.0);
				else fprintf(fforet_avg, "%lf\t%lf\t%lf\n",tts[t],etas_mu*(tts[t]-tts[t-1])+nev[t]/(1.0*NgridT),etas_mu+rev[t]/(1.0*NgridT));
			}
			fclose(fforet_avg);
		}

		if (print_forex0) {
//			convert_geometry(crst, ev_x, &ev_x_new, 1, 0);
			for (int n=1; n<=NgridT_out; n++) {
				ev_x[n]+=etas_mu*(tts[Ntts]-tts[0]);
				ev_x[n]*=1.0/NgridT;
			}
			csep_forecast(print_forex, crst, ev_x, 1);
		}

		if (print_LL){
			tt0=tts[0];
			tt1=tts[Ntts];
			current_main=0;
			Ldum=0.0;
			tnow=tt0;
			int j0=1;
			while (current_main<Nm && eqkfm[current_main].t<tt0) current_main++;
			while (current_main<Nm && eqkfm[current_main].t<tt1){
				if (tnow<eqkfm[current_main].t){
					for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<eqkfm[current_main].t) {
						fprintf(fLLev,"%.10e \t%.5e \t%.5e \t%.5e \t %.5e\n",cat.t[j], cat.lat0[j], cat.lon0[j], cat.depths0[j], log(etas_mu+rate[j]));
					}
					j0+=N;
				}
				tnow=eqkfm[current_main].t+tw;
				current_main+=1;
			}
			if (tnow<tt1){
				for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<tt1) {
					fprintf(fLLev,"%.10e \t%.5e \t%.5e \t%.5e \t %.5e\n",cat.t[j], cat.lat0[j], cat.lon0[j], cat.depths0[j], log(etas_mu+rate[j]));
				}
			}

			fclose(fLLev);
		}
	}

	free_dvector(dumrate,1,cat.Z);
	free_dvector(rate, 1,cat.Z);
	free_dvector(ev_x,1,NgridT);
	free_dvector(nev,1,Ntts);
	free_dvector(rev,1,Ntts);
	return(err);

}


int forecastETAS(struct catalog cat, struct pscmp *DCFS, double tt0, double tt1, int points[], double *out_NeX, double *NeT,
		double *Rate_end, int N, int Neqks, double *R, double etas_c, double etas_p, double etas_K){
//assumes stress grows linearly during each time step.
//NB assumes that times, cat, DCFS, NTS are always the same (each time function is called).
// if points==NULL, use sequence [1,2,3,...,N].
// cmbdata can be NULL, and will be ignored.
// if backrate==1, will assume it's 1 for all points. If not, it should be an array containing NgridT times the ratio between avg rate and rate at that point. (so that the sum of back_rate is NgridT).
	//fixme check: indices of R[0...cat.Z-1] or {[1...cat.Z]

  double t, t0;
  int n;
  static double **events;	//will contain times, magnitudes of all earthquakes (both from catalog and from DCFS).
  int Neq, err=0, errtot=0;
  int is_incat, is_inDCFS;
  int next_eqk, counter_eqk;
  int cat_i, DCFS_i;
  double *NeX;
  double **Rprivate;
  static double *NeXdum, *ReX;
  static int *TS_eqk;
  static int firsttimein=1;
  static int **which_eqk;	//contains list of eqk for each grid point, and index of the point within DCFS;
  static int *num_eqk;
  static float ** cat_weights;
  static int ** DCFS_whichpt;
  int reach_end;
  int nthreads, nthreadstot=omp_get_max_threads();
  static int **indices;	//contains indices of events referred to DCFS and cat(with offset of 1 - see below).

  if (firsttimein==1){
		firsttimein=0;
		//todo: use previous information to match DCFS and cat.
		events=union_cats(cat.t+1, timesfrompscmp(DCFS,Neqks), cat.mag+1, magsfrompscmp(DCFS,Neqks), cat.Z, Neqks, 1e-6, 1.0, &indices, &Neq);	//todo parameters should not be hardwired.
		NeXdum=dvector(1,N);
		ReX=dvector(1,N);

		//list of all events, total no. of events affecting each point (Neq: from combining DCFS and cat).
		which_eqk=imatrix(1,N,0,Neq-1);
		num_eqk=ivector(1,N);
		cat_weights=matrix(1,N,0,Neq-1);
		DCFS_whichpt=imatrix(1,N,0,Neq-1);
		TS_eqk=ivector(0,Neq-1);

		for (int n=1; n<=N; n++) num_eqk[n]=0;

		for(int eq=0;eq<Neq;eq++){
			int counter1, counter2;
			counter1=counter2=1;
			cat_i=indices[1][eq]+1;	//	NB:cat_i[i]==0 means no events selected.
			DCFS_i=indices[2][eq];  //	NB:DCFS_i[i]==-1 means no events selected.

			if(DCFS_i>0){
				if (DCFS[DCFS_i].index_cat!=cat_i) {
					printf("** Error!! mismatch in DCFS and catalog indices!!** \n");
					if (flog){
						fprintf(flog, "** Error!! mismatch in DCFS and catalog indices!! (t=%.3lf\tmag=%.3lf)** \n", DCFS[DCFS_i].t, DCFS[DCFS_i].m);
					}
				}
			}

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
  }

  for (int z=1; z<=cat.Z; z++) if (cat.t[z]>=tt0 && cat.t[z]<tt1) R[z]=0.0;
  NeX= (out_NeX)? out_NeX : NeXdum;

  if (Rate_end) for(int m=1;m<=N;m++) ReX[m]=0.0;
  if (NeX) for(int m=1;m<=N;m++) NeX[m]=0.0;

  if (NeT) *NeT=0;
  if (Rate_end) *Rate_end=0;

//  printf ("max no of threads=%d.\n", omp_get_max_threads());

  Rprivate=dmatrix(0,nthreadstot-1, 0, cat.Z);
  for (int t=0; t<omp_get_max_threads(); t++){
	  for (int eq=0; eq<=cat.Z; eq++) Rprivate[t][eq]=0.0;
  }

//#pragma omp parallel for private(err, n, next_eqk, reach_end, cat_i, DCFS_i, counter_eqk) reduction(+:errtot)
  for(int m=1;m<=N;m++){

	nthreads=omp_get_num_threads();

	n=(points==0)? m : points[m];
    counter_eqk=0;
    reach_end=0;

    err=0;
    ReX[m]=0.0;

    while (counter_eqk<num_eqk[n] && events[1][which_eqk[n][counter_eqk]]<tt0) counter_eqk+=1;
	next_eqk=which_eqk[n][counter_eqk];

	while (counter_eqk<num_eqk[n] && events[1][next_eqk]<tt1){
		cat_i= indices[1][next_eqk]+1;
			if (cat_i!=0) {
			for (int eq_prec=0; eq_prec<counter_eqk; eq_prec++){
				next_eqk=which_eqk[n][eq_prec];
				DCFS_i= indices[2][next_eqk];
				t=cat.t[cat_i]-DCFS[DCFS_i].t;
				Rprivate[omp_get_thread_num()][cat_i]+= cat_weights[n][counter_eqk]*etas_K*DCFS[eq_prec].cmb[DCFS_whichpt[n][counter_eqk]]*pow((t+etas_c),-1.0*etas_p);
			}
		}
		counter_eqk+=1;
		next_eqk=which_eqk[n][counter_eqk];
	}

	for (int eq_prec=0; eq_prec<counter_eqk; eq_prec++){
		next_eqk=which_eqk[n][eq_prec];
		DCFS_i= indices[2][next_eqk];
		t=tt1-DCFS[DCFS_i].t;
		t0=fmax(0.0, tt0-DCFS[DCFS_i].t);
		ReX[m]+= etas_K*DCFS[eq_prec].cmb[DCFS_whichpt[n][counter_eqk]]*pow((t+etas_c),-1.0*etas_p);
		if (etas_p==1.0) NeX[m]+= etas_K*DCFS[next_eqk].cmb[DCFS_whichpt[n][counter_eqk]]*(log(t+etas_c)-log(t0+etas_c));
		else NeX[m]+= (1.0/(1.0-etas_p))*etas_K*DCFS[next_eqk].cmb[DCFS_whichpt[n][counter_eqk]]*(pow((t+etas_c),1.0-etas_p)-pow((t0+etas_c),1.0-etas_p));
	}
}


  for (int t=0; t<nthreads; t++){
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
