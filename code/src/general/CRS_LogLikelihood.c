/*
 * calculate_LogLikelihood.c
 *
 *  Created on: Jul 22, 2013
 *      Author: camcat
 */

#include "CRS_LogLikelihood.h"

int CRSforecast (double *LL, int Nsur, int Nslipmod, struct pscmp *DCFS, struct eqkfm *eqkfm_aft, struct eqkfm *eqkfm0, struct eqkfm *eqkfm1, struct flags flags,
		double *tevol, struct crust crst, struct Coeff_LinkList *AllCoeff, int NTScont, int NTSdisc, int Nm, int NgridT, double **focmec, int *fmzonelim, int NFM,
		long *seed, struct catalog cat, double *times, double tstart, double *tts, int Ntts, double tw, double Asig, double ta, double r0,
		double **all_gammas0, int multiple_input_gammas, int fromstart, double Hurst,
		char * print_cmb, char *print_forex, char *print_foret, char * printall_cmb, char *printall_forex, char *printall_foret, char *print_LL){


	//Similar to CRSLogLikelihood, but loops over time steps to produce time forecast.
	//recfault= [0,1,2] means: don't vary rec. fault, vary (choose random one), vary and sample all catalog in order.
	//tstart=starting time of entire period considered (may be <tt0). only used if fromstart==1;
	//if multiple_input_gammas==1, expect allgammas0[1...Nsur][1...NgridT]; else, allgammas0is pointer to vector: [1...NgridT].
	//same for multiple_input_gammas (if multiple_output_gammas, return gammas that would give average rate).

	/*INPUT:-------------------------------------------------------------//
	 *
	 * Nsur= no. of iterations
	 * DCFS= array containing stress steps.
	 * eqkfm_aft= array of afterslip:
	 * 		if splines==1
	 * 		if splines==0
	 * eqkfm0= array of mainshocks;
	 * eqkfm1= array of aftershocks (sources);
	 * tevol= if splines==0, contains factor by which afterslip snapshot is rescaled at each time.
	 * crst= general model info.
	 * AllCoeff= Okada Coefficients for all mainshock slip models;
	 * NTScont=no. of continuous time steps (size of tevol, times2)
	 * NTSdisc=no of snapshots (size of DCFS);  size of eqkfm1
	 * Nm=no. of mainshocks;	size of eqkfm0
	 * focmec= array of focal mechanisms parameters.
	 * NFM= no. of focal mechanisms
	 *
	 * flags:
	 * afterslip, aftershocks;
	 * vary_recfault: 0=use fix planes, 1=use foc. mec, 2=ose OOPs;
	 * vary_slipmodel
	 * gridpoints_err
	 * splines
	 *
	 * struct catalog cat= catalog for which LL ir calculated.
	 * times: time steps for numerical integration.
	 * tstart: overall time at which calculation of rates should start;
	 * tts: list of time steps (tts[0] is starting time; tts[1...Ntts] are times at which rate is calculated).
	 * tw: time window to be ignored after each mainshock;
	 *
	 * RS parameters and flags:
	 * Asig, ta, r0;
	 * all_gammas0=initial gammas;	is NULL, use gamma=ta/Asig;
	 *
	 * multiple_input_gammas: indicates if a set of gammas per iteration is given;
	 * multiple_output_gammas: indicates if a set of gammas per iteration should be provided as output;
	 * fromstart: indicates if the gammas refer to tstart (fromstart=1) or t0(fromstart=0).
	 *
	 * the following are filenames to which output is written (ignored if NULL is passed). all with extension, except printall_foret(will add _cumu.dat, _rate.dat)
	 * printall_cmb, printall_forex: prints out a file containing the values at each gridpoint for all iterations (cmb value is for DCFS[0]).
	 * printall_foret: prints out a file containing time evolution of forecast for all iterations.
	 * print_LL: prints out a file containing log(r_ek) for all events ek.
	 *
	 */

	double Ldum_out, Nev;

	char fname[120];
	static double **DCFSrand;
	static double *dumrate, *gammas, *rates_x, *rate, *ev_x, *ev_x_new=0;
	double sum, sum1, sum2, integral;
	double Ldum;
	double fin_rate;
	double *gammas0;
	double *nev, *rev, *nev_avg, *rev_avg, *ev_x_avg, *ev_x_dum;
	double *cmb, *cmb_avg;
	int N, NgridT_out= crst.uniform ? (crst.nLat_out*crst.nLon_out*crst.nD_out) : NgridT;
	int err, Nsur_over_Nslipmod=Nsur/Nslipmod;
	int uniform_bg_rate=0;
	int current_main;
	double tnow, tt0, tt1;
	FILE *fforex, *fcmb, *fforet1, *fforet2, *fLLev;
	FILE *fforex_avg, *fcmb_avg, *fforet_avg;

	if (all_gammas0==NULL)	uniform_bg_rate=1;


	DCFSrand= (flags.afterslip) ? dmatrix(0,NTScont,1,NgridT) : NULL;
	dumrate=dvector(1,cat.Z);
	rate=dvector(1,cat.Z);
	gammas=dvector(1,NgridT);
	rates_x=dvector(1,NgridT);
	ev_x=dvector(1,NgridT);
	ev_x_avg=dvector(1,NgridT);
	ev_x_dum=dvector(1,NgridT);
	cmb=dvector(1,NgridT);
	cmb_avg=dvector(1,NgridT);
	nev=dvector(1,Ntts);
	rev=dvector(1,Ntts);
	nev_avg=dvector(1,Ntts);
	rev_avg=dvector(1,Ntts);
	for (int n=1; n<=NgridT; n++) {
		ev_x_avg[n]=0.0;
		cmb_avg[n]=0.0;
	}
	for (int t=1; t<=Ntts; t++) {
		nev_avg[t]=0.0;
		rev_avg[t]=0.0;
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
	if (print_LL) fLLev=fopen(print_LL,"w");
	if (print_cmb) fcmb_avg=fopen(print_cmb,"w");
	if (print_forex) fforex_avg=fopen(print_forex,"w");
	if (print_foret) fforet_avg=fopen(print_foret,"w");
	if (printall_cmb) fcmb=fopen(printall_cmb,"w");
	if (printall_forex) fforex=fopen(printall_forex,"w");
	if (printall_foret) {
		sprintf(fname, "%s_cumu.dat",printall_foret);
		fforet1=fopen(fname,"w");
		sprintf(fname, "%s_rates.dat",printall_foret);
		fforet2=fopen(fname,"w");
		//for (int t=1; t<=Ntts; t++) fprintf(fforet, "%lf\t",tts[t]);
		//fprintf(fforet, "\n");
	}
	for (int nsur=1; nsur<=Nsur; nsur++){

		for (int n=1; n<=NgridT; n++) ev_x[n]=0.0;
		for(int i=1;i<=cat.Z;i++) dumrate[i]=0.0;

		printf("%d...",nsur);
		fflush(stdout);
		if (all_gammas0) gammas0= (multiple_input_gammas)? all_gammas0[nsur] : *all_gammas0;
		flags.which_recfault= (NFM!=0 && Nsur%NFM==0)? nsur%NFM+1 : 0;	//which_recfault=0 means: choose random one.
		flags.new_slipmodel= !(nsur % Nsur_over_Nslipmod);

		tt0=tts[0];
		tt1=tts[Ntts];
		//Set starting rates:
		if (fromstart){
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, eqkfm1, flags, tevol, times, Nm, crst, AllCoeff, NTScont, NTSdisc, focmec, fmzonelim, NFM, seed, (int *) 0, tstart, tt1, Hurst, nsur==1);
			for (int n=1; n<=NgridT; n++) gammas[n]= (uniform_bg_rate)? ta/Asig : gammas0[n];
			err=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tstart, tt0, Asig, ta, (int *) 0, (double *) 0, (double *) 0, (double *) 0, NgridT,
					NTScont, NTSdisc, gammas,	dumrate, 1);
		}
		else{
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, eqkfm1, flags, tevol, times, Nm, crst, AllCoeff, NTScont, NTSdisc, focmec, fmzonelim, NFM, seed, (int *) 0, tt0, tt1, Hurst, nsur==1);
			for (int n=1; n<=NgridT; n++) gammas[n]= (uniform_bg_rate)? ta/Asig : gammas0[n];
		}

		current_main=0;
		tnow=tts[0];
			for (int t=1; t<=Ntts; t++){
			//Calculate seismicity evolution (skipping a time window after each mainshock):
			tt0=tts[t-1];
			tt1=tts[t];
			while (current_main<Nm && eqkfm0[current_main].t<tt0) current_main++;
			while (current_main<Nm && eqkfm0[current_main].t<tt1){
				if (tnow<eqkfm0[current_main].t){
					err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tnow, eqkfm0[current_main].t, Asig, ta, 0, ev_x_dum, &sum, &fin_rate, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
					for (int i=1; i<=NgridT; i++) ev_x[i]+=ev_x_dum[i];
					nev[t]=sum;
					rev[t]=fin_rate;
					nev_avg[t]+=(sum)/(1.0*Nsur);
					rev_avg[t]+=(fin_rate)/(1.0*Nsur);
					err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, eqkfm0[current_main].t, eqkfm0[current_main].t+tw, Asig, ta, 0, 0, &sum, 0, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
					tnow=eqkfm0[current_main].t+tw;
				}
				else if (tnow<eqkfm0[current_main].t+tw){
					err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tnow, eqkfm0[current_main].t+tw, Asig, ta, 0, 0, &sum, 0, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
					tnow=eqkfm0[current_main].t+tw;
				}
				current_main+=1;
			}
			if (tnow<tt1){
				err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tnow, tt1, Asig, ta, 0, ev_x_dum, &sum, &fin_rate, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
				for (int i=1; i<=NgridT; i++) ev_x[i]+=ev_x_dum[i];
				nev[t]=sum;
				rev[t]=fin_rate;
				nev_avg[t]+=(sum)/(1.0*Nsur);
				rev_avg[t]+=(fin_rate)/(1.0*Nsur);
				tnow=tt1;
			}
			for(int i=1;i<=cat.Z;i++) if(cat.t[i]>tt0 && cat.t[i]<=tt1) rate[i]+=dumrate[i]/(1.0*Nsur);
		}
		if (err==1) break;

		for (int i=1; i<=NgridT; i++) ev_x_avg[i]+=ev_x[i]/(1.0*Nsur);

		if (print_cmb || printall_cmb){
			sum_DCFS(DCFS, &cmb, NTSdisc, NgridT);
			if (print_cmb) for (int i=1; i<=NgridT; i++) cmb_avg[i]+=cmb[i]*(1.0/Nsur);
			if (printall_cmb){
				convert_geometry(crst,cmb, &ev_x_new, 0, 0);
				for (int n=1; n<=NgridT_out; n++) fprintf(fcmb, "%lf\t", ev_x_new[n]);
				if (nsur <Nsur) fprintf(fcmb, "\n");
				fflush(fcmb);
			}
		}

		if (printall_forex) {
			convert_geometry(crst, ev_x, &ev_x_new, 1, 0);
			for (int n=1; n<=NgridT_out; n++) fprintf(fforex, "%lf\t", ev_x_new[n]*r0/NgridT);
			if (nsur <Nsur) fprintf(fforex, "\n");
			fflush(fforex);
		}

		if (printall_foret) {
			for (int t=1; t<=Ntts; t++) {
				fprintf(fforet1, "%lf\t",nev[t]*r0/NgridT);
				fprintf(fforet2, "%lf\t",rev[t]*r0/NgridT);
			}
			if (nsur <Nsur) fprintf(fforet1, "\n");
			if (nsur <Nsur) fprintf(fforet2, "\n");
			fflush(fforet1);
			fflush(fforet2);
		}

	}

	//calculate average rate and LL:
	if (!err){

		if (printall_cmb) fclose(fcmb);
		if (printall_forex) fclose(fforex);
		if (printall_foret) {
			fclose(fforet1);
			fclose(fforet2);
		}

		if (print_foret) {
			for (int t=1; t<=Ntts; t++) fprintf(fforet_avg, "%lf\t%lf\t%lf\n",tts[t],nev_avg[t]*r0/(1.0*NgridT),rev_avg[t]*r0/(1.0*NgridT));
			fclose(fforet_avg);
		}

		if (print_forex) {
			convert_geometry(crst, ev_x_avg, &ev_x_new, 1, 0);
			for (int n=1; n<=NgridT_out; n++) ev_x_new[n]*=r0/NgridT;
			csep_forecast(print_forex, crst, ev_x_new);
		}

		if (print_cmb) {
			convert_geometry(crst, cmb_avg, &ev_x_new, 0, 0);
			csep_forecast(print_cmb, crst, ev_x_new);
		}
		if (print_LL || LL){
			tt0=tts[0];
			tt1=tts[Ntts];
			current_main=0;
			Ldum=0.0;
			tnow=tt0;
			int j0=1;
			while (current_main<Nm && eqkfm0[current_main].t<tt0) current_main++;
			while (current_main<Nm && eqkfm0[current_main].t<tt1){
				if (tnow<eqkfm0[current_main].t){
					for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<eqkfm0[current_main].t) {
						Ldum+=log(r0*rate[j]);
						if (print_LL) fprintf(fLLev,"%.10e \t %.5e\n",cat.t[j], log(r0*rate[j]));
					}
					j0+=N;
				}
				tnow=eqkfm0[current_main].t+tw;
				current_main+=1;
			}
			if (tnow<tt1){
				for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<tt1) {
					Ldum+=log(r0*rate[j]);
					if (print_LL) fprintf(fLLev,"%.10e \t %.5e\n",cat.t[j], log(r0*rate[j]));
				}
			}

			if (print_LL) fclose(fLLev);
			if (LL){
				integral=0.0;
				for (int t=1; t<=Ntts; t++) integral+= nev_avg[t];
				*LL=Ldum-integral*r0/(1.0*NgridT);
			}
		}
	}

	if (flags.afterslip) free_dmatrix(DCFSrand, 0,NTScont,1,NgridT);
	free_dvector(dumrate,1,cat.Z);
	free_dvector(rate, 1,cat.Z);
	free_dvector(gammas, 1,NgridT);
	free_dvector(rates_x, 1,NgridT);
	free_dvector(ev_x,1,NgridT);
	free_dvector(ev_x_avg,1,NgridT);
	free_dvector(cmb,1,NgridT);
	free_dvector(cmb_avg,1,NgridT);
	free_dvector(nev,1,Ntts);
	free_dvector(rev,1,Ntts);
	free_dvector(nev_avg,1,Ntts);
	free_dvector(rev_avg,1,Ntts);
	return(err);

}

int CRSLogLikelihood (double *LL, double *Ldum0_out, double *Nev, double *I, double *r_out, int Nsur, int Nslipmod, struct pscmp *DCFS,
		struct eqkfm *eqkfm_aft, struct eqkfm *eqkfm0, struct eqkfm *eqkfm1, struct flags flags, double Hurst,
		double *tevol, struct crust crst, struct Coeff_LinkList *AllCoeff, int NTScont, int NTSdisc, int Nm, int NgridT, double **focmec, int *fmzonelim, int NFM,
		long *seed, struct catalog cat, double *times, double tstart, double tt0, double tt1, double tw, double Asig, double ta, double r0, int fixr,
		double **all_gammas0, double **all_new_gammas, int multiple_input_gammas, int multiple_output_gammas, int fromstart,
		char * printall_cmb, char *printall_forex, int refresh){

	//recfault= [0,1,2] means: don't vary rec. fault, vary (choose random one), vary and sample all catalog in order.
	//tstart=starting time of entire period considered (may be <tt0). only used if fromstart==1;
	//if multiple_input_gammas==1, expect allgammas0[1...Nsur][1...NgridT]; else, allgammas0is pointer to vector: [1...NgridT].
	//same for multiple_input_gammas (if multiple_output_gammas, return gammas that would give average rate).

	/*INPUT:-------------------------------------------------------------//
	 *
	 * Nsur= no. of iterations
	 * DCFS= array containing stress steps.
	 * eqkfm_aft= array of afterslip:
	 * 		if splines==1
	 * 		if splines==0
	 * eqkfm0= array of mainshocks;
	 * eqkfm1= array of aftershocks (sources);
	 * tevol= if splines==0, contains factor by which afterslip snapshot is rescaled at each time.
	 * crst= general model info.
	 * AllCoeff= Okada Coefficients for all mainshock slip models;
	 * NTScont=no. of continuous time steps (size of tevol, times2)
	 * NTSdisc=no of snapshots (size of DCFS);  size of eqkfm1
	 * Nm=no. of mainshocks;	size of eqkfm0
	 * focmec= array of focal mechanisms parameters.
	 * NFM= no. of focal mechanisms
	 *
	 * flags:
	 * afterslip, aftershocks;
	 * vary_recfault: 0=use fix planes, 1=use foc. mec, 2=ose OOPs;
	 * vary_slipmodel
	 * gridpoints_err
	 * splines
	 *
	 * struct catalog cat= catalog for which LL ir calculated.
	 * times: time steps for numerical integration.
	 * tstart: overall time at which calculation of rates should start;
	 * [tt0, tt1]: time for which LL is calculated;
	 * tw: time window to be ignored after each mainshock;
	 *
	 * RS parameters and flags:
	 * Asig, ta, r0, fixr;
	 * all_gammas0=initial gammas;	is NULL, use gamma=ta/Asig;
	 *
	 * multiple_input_gammas: indicates if a set of gammas per iteration is given;
	 * multiple_output_gammas: indicates if a set of gammas per iteration should be provided as output;
	 * fromstart: indicates if the gammas refer to tstart (fromstart=1) or t0(fromstart=0).
	 *
	 * the following are filenames to which output is written (ignored if NULL is passed).
	 * printall_cmb, printall_forex: prints out a file containing the values at each gridpoint for all iterations (cmb value is for DCFS[0]).
	 *
	// OUTPUT:-------------------------------------------------------------//
	 *
	 * Output variables are incremented by the amount calculated for relevant time period; if fixrate=0, LL is recalculated from the start.
	 *
 	 * LL= LogLikelihood;
 	 * LL= LLdum0+Nlog(r)-r*I;
	 * Ldum0= summation part of LogLikelihood (w/o rate dependence);
	 * I= tot no. events (integral part of LL) / r.;
	 * all_new_gammas= final gammas.
	 *
	 */

	static double **DCFSrand;
	static double *dumrate, *gammas, *rates_x, *rate, *ev_x;
	double sum, sum1, sum2, integral, Ldum0, r;
	double *gammas0;
	static int first_timein=1, N;
	double Ntot, Itot, LLdum0tot;
	int err, Nsur_over_Nslipmod=Nsur/Nslipmod;
	int uniform_bg_rate=0;
	int current_main, j0;
	double tnow;
	FILE *fforex, *fcmb;

	if (LL && verbose_level>0) printf("Calculating LL for Asig=%lf, ta=%lf ...", Asig, ta);

	if (!all_gammas0)	uniform_bg_rate=1;

	if (first_timein==1){
		if (flog) fprintf(flog,"Setting up variables in CRSLogLikelihood...\n"); fflush(flog);
		first_timein=0;
		DCFSrand= (flags.afterslip) ? dmatrix(0,NTScont,1,NgridT) : NULL;
		dumrate=dvector(1,cat.Z);
		rate=dvector(1,cat.Z);
		gammas=dvector(1,NgridT);
		rates_x=dvector(1,NgridT);
		ev_x=dvector(1,NgridT);
	}

	sum=sum1=sum2=0.0;
	integral=0.0;
	for(int i=1;i<=cat.Z;i++) rate[i]=0.0;


	if (all_new_gammas!=0 && (multiple_output_gammas==0)) for (int n=1; n<=NgridT; n++) rates_x[n]=0.0;	//need to save average rate at each location, to find final gamma giving average rate.
	//for (int ndt=1; ndt<=NDT; ndt++) net[ndt]=0.0;

	sum=sum1=sum2=0;
	err=0;
	if (printall_cmb) fcmb=fopen(printall_cmb,"w");
	if (printall_forex) fforex=fopen(printall_forex,"w");
	for (int nsur=1; nsur<=Nsur; nsur++){

		//for(int i=1;i<=cat.Z;i++) dumrate[i]=0.0;

		if (all_gammas0) gammas0= (multiple_input_gammas)? all_gammas0[nsur] : *all_gammas0;
		flags.which_recfault= (NFM!=0 && Nsur%NFM==0)? nsur%NFM+1 : 0;	//which_recfault=0 means: choose random one.
		flags.new_slipmodel= !(nsur % Nsur_over_Nslipmod);

		//Set starting rates:
		if (fromstart){
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, eqkfm1, flags, tevol, times, Nm, crst, AllCoeff, NTScont, NTSdisc, focmec, fmzonelim, NFM, seed, (int *) 0, tstart, tt1, Hurst, refresh && nsur==1);
			for (int n=1; n<=NgridT; n++) gammas[n]= (uniform_bg_rate)? ta/Asig : gammas0[n];
			err=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tstart, tt0, Asig, ta, (int *) 0, (double *) 0, (double *) 0, (double *) 0, NgridT,
					NTScont, NTSdisc, gammas,	dumrate, 1);
		}
		else{
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, eqkfm1, flags, tevol, times, Nm, crst, AllCoeff, NTScont, NTSdisc, focmec, fmzonelim, NFM, seed, (int *) 0, tt0, tt1, Hurst, refresh && nsur==1);
			for (int n=1; n<=NgridT; n++) gammas[n]= (uniform_bg_rate)? ta/Asig : gammas0[n];
		}

		//Calculate seismicity evolution (skipping a time window after each mainshock):
		current_main=0;
		tnow=tt0;
		while (current_main<Nm && eqkfm0[current_main].t<tt0) current_main++;
		while (current_main<Nm && eqkfm0[current_main].t<tt1){
			if (tnow<eqkfm0[current_main].t){
				err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tnow, eqkfm0[current_main].t, Asig, ta, 0, 0, &sum, 0, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
				integral+=(sum)/(1.0*Nsur);
				err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, eqkfm0[current_main].t, eqkfm0[current_main].t+tw, Asig, ta, 0, 0, &sum, 0, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
				tnow=eqkfm0[current_main].t+tw;
			}
			else if (tnow<eqkfm0[current_main].t+tw){
				err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tnow, eqkfm0[current_main].t+tw, Asig, ta, 0, ev_x, &sum, 0, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
				tnow=eqkfm0[current_main].t+tw;
			}
			current_main+=1;
		}
		if (tnow<tt1){
			err+=forecast_stepG2_new(cat, times, DCFSrand, DCFS, tnow, tt1, Asig, ta, 0, ev_x, &sum, 0, NgridT, NTScont, NTSdisc, gammas, dumrate, 1);
			integral+=(sum)/(1.0*Nsur);
		}

		if (err==1) break;

		for(int i=1;i<=cat.Z;i++) if(cat.t[i]>=tt0 && cat.t[i]<tt1) rate[i]+=1.0*dumrate[i]/(1.0*Nsur);
		if (all_new_gammas) {
			if (multiple_output_gammas) for (int n=1; n<=NgridT; n++) all_new_gammas[nsur][n]=gammas[n];
			else for (int n=1; n<=NgridT; n++) rates_x[n]+=(1.0/gammas[n])/(1.0*Nsur);
		}

		if (printall_cmb) {
			for (int n=1; n<=NgridT; n++) fprintf(fcmb, "%lf\t", DCFS[0].cmb[n]);
			if (nsur <Nsur) fprintf(fcmb, "\n");
		}
		if (printall_forex) {
			for (int n=1; n<=NgridT; n++) fprintf(fforex, "%lf\t", ev_x[n]);
			if (nsur <Nsur) fprintf(fforex, "\n");
		}
	}

	//calculate average rate and LL:
	if (!err){
		if (printall_cmb) fclose(fcmb);
		if (printall_forex) fclose(fforex);

		Ldum0=0.0;
		N=0;
		current_main=0;
		tnow=tt0;
		j0=1;
		while (current_main<Nm && eqkfm0[current_main].t<tt0) current_main++;
		while (current_main<Nm && eqkfm0[current_main].t<tt1){
			if (tnow<eqkfm0[current_main].t){
				for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<eqkfm0[current_main].t) {
					N+=1;
					Ldum0+=log(rate[j]);
				}
				j0+=N;
			}
			tnow=eqkfm0[current_main].t+tw;
			current_main+=1;
		}
		if (tnow<tt1){
			for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<tt1) {
				N+=1;
				Ldum0+=log(rate[j]);
			}
		}

		Ntot= (Nev) ? N+*Nev : N;
		Itot= (I)? integral/NgridT + *I : integral/NgridT;
		LLdum0tot= (Ldum0_out) ? Ldum0+*Ldum0_out : Ldum0;
		r = (fixr)? r0 : Ntot/Itot;
		if(r==0.0 && verbose_level>0) printf("ERROR: ta=%lf  Asig=%lf r=%e\n",ta,Asig,r);

		if (LL) *LL=LLdum0tot+ Ntot*log(r) - r*Itot;
		if (Ldum0_out) *Ldum0_out=LLdum0tot;
		if (I) *I=Itot;
		if (Nev) *Nev=Ntot;
		if (r_out) *r_out=r;
		if (all_new_gammas!=0 && multiple_output_gammas==0) for (int n=1; n<=NgridT; n++) (*all_new_gammas)[n]=1.0/rates_x[n];

		if (LL && verbose_level>0) printf("LL=%lf\n", *LL);
	}

	if (flog && first_timein) printf(flog,"done.\n");
	return(err);

}
