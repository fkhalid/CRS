/*
 * calculate_LogLikelihood.c
 *
 *  Created on: Jul 22, 2013
 *      Author: camcat
 */

#include "CRS_LogLikelihood.h"

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../geom/convert_geometry.h"
#include "../inp_out/print_output.h"
#include "../inp_out/write_csep_forecast.h"
#include "../util/error.h"
#include "../util/nrutil.h"
#include "calculateDCFSperturbed.h"
#include "forecast_stepG.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int CRSforecast(double *LL, int Nsur, struct pscmp *DCFS, struct eqkfm *eqkfm_aft,
				struct eqkfm *eqkfm0, struct flags flags,
				struct crust crst, struct Coeff_LinkList *AllCoeff,
				int NTScont, int Nm, int Na, int NgridT, double **focmec, int *fmzonelim,
				int NFM, long *seed, struct catalog cat, double *times, double tstart, double tt0, double tt1,
				double dtstep, double Asig, double ta, double r0, double **all_gammas0,
				int multiple_input_gammas, int fromstart, char * print_cmb0,
				char *print_forex0, char *print_foret, char * printall_cmb, char *printall_forex,
				char *printall_foret, char *print_LL, int refresh) {

	//Similar to CRSLogLikelihood, but loops over time steps to produce time forecast.
	//recfault= [0,1,2] means: don't vary rec. fault, vary (choose random one), vary and sample all catalog in order.
	//tstart=starting time of entire period considered (may be <tt0). only used if fromstart==1;
	//if multiple_input_gammas==1, expect allgammas0[1...Nsur][1...NgridT]; else, allgammas0is pointer to vector: [1...NgridT].
	// this can be useful if non stationary values of gamma are used (...)
	//same for multiple_input_gammas (if multiple_output_gammas, return gammas that would give average rate).
	//all_gammas0[0] can be an array of starting values for gamma; they will be non steady-state values. Otherwise, ...

	/*INPUT:-------------------------------------------------------------//
	 *
	 * Nsur= no. of iterations
	 * DCFS= array containing stress steps.
	 * eqkfm_aft= array of afterslip:
	 * 		if splines==1
	 * 		if splines==0
	 * eqkfm0= array of mainshocks;
	 * crst= general model info.
	 * AllCoeff= Okada Coefficients for all mainshock slip models;
	 * NTScont=no. of continuous time steps (size of times2)
	 * Nm=no. of mainshocks;	size of eqkfm0
	 * Na= no. of events with afterslip
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

	//fixme: v important. if LLinversion=0 and multiple slip models are given, something equivalent to "refresh" should be done here (see CRSLoglikelihood).

	// [Fahad] Variables used for MPI.
	int procId = 0, numProcs = 1;
	int start, end, partitionSize;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	#endif

	double Ldum_out, Nev;

	char fname[120];
	char print_forex[120], print_cmb[120], print_cmbpost[120];
	static double **DCFSrand;
	static double *dumrate, *gammas, *rate, *ev_x, *ev_x_new=NULL, *flat_grid;
	double integral;
	double Ldum;
	double *gammas0;
	double *nev, *rev, *nev_avg, *rev_avg, *ev_x_avg, *ev_x_pre, *ev_x_dum;
	double **nev_allsnapshots=NULL, **nev_allsnapdum=NULL;
	double *cmb, *cmb_avg, *cmbpost, *cmbpost_avg;
	int N, NgridT_out= crst.uniform ? (crst.nLat_out*crst.nLon_out*crst.nD_out) : NgridT;
	int err;
	int msnap;
	int uniform_bg_rate=0;
	int current_main, which_recfault;
	int NgridTsnaps;
	int flatten_allsnaps=1;	//boolean to decide whether all snapshots should be flattened (i.e. sum over depth layers) to avoid huge files.
	int Ntts=ceil((tt1-tt0)/dtstep);
	double tnow;
	FILE *fforex, *fcmb, *fforet1, *fforet2, *fLLev;
	FILE *fforet_avg, *foutallsnap;

	int all_snapshots=1;	//todo pass to function.

	if (all_gammas0==NULL)	uniform_bg_rate=1;

	//Initialize stress fields:
	DCFSrand= (flags.afterslip) ? dmatrix(0,NTScont,1,NgridT) : NULL;
	if (flags.afterslip){
		DCFSrand= dmatrix(0,NTScont,1,NgridT);
		for (int i=0; i<NTScont; i++){
			for (int j=1; j<=NgridT; j++){
				DCFSrand[i][j]=0;
			}
		}
	}

	else DCFSrand=NULL;

	dumrate=dvector(1,cat.Z);
	rate=dvector(1,cat.Z);
	gammas=dvector(1,NgridT);
	ev_x=dvector(1,NgridT);
	ev_x_avg=dvector(1,NgridT);
 	ev_x_pre=dvector(1,NgridT);
	ev_x_dum=dvector(1,NgridT);
	ev_x_new=dvector(1,NgridT_out);
	cmb=dvector(1,NgridT);
	cmbpost=dvector(1,NgridT);
	cmbpost_avg=dvector(1,NgridT);
	cmb_avg=dvector(1,NgridT);
	nev=dvector(1,Ntts);
	rev=dvector(1,Ntts);
	nev_avg=dvector(1,Ntts);
	rev_avg=dvector(1,Ntts);

	if (all_snapshots){
		if (flatten_allsnaps && !crst.uniform){
			print_screen("Warning: can not flatten snapshots over depth, since grid is not uniform.\n");
			print_logfile("Warning: can not flatten snapshots over depth, since grid is not uniform.\n");
			flatten_allsnaps=0;
		}
		NgridTsnaps=(flatten_allsnaps)? NgridT_out/crst.nD_out : NgridT_out;
		flat_grid=dvector(1, NgridTsnaps);
		nev_allsnapshots= dmatrix(1,Ntts, 1, NgridT);
		nev_allsnapdum=dmatrix(0,Ntts-1, 1, NgridT);	//indices are like this because of how array is address in rate_state_evolution.
	}

	for (int n=1; n<=NgridT; n++) {
		ev_x_avg[n]=0.0;
		ev_x_pre[n]=0.0;
		cmb_avg[n]=0.0;
		cmbpost_avg[n]=0.0;
	}
	for (int t=1; t<=Ntts; t++) {
		nev_avg[t]=0.0;
		rev_avg[t]=0.0;
		nev[t]=0.0;
		rev[t]=0.0;
	}
	N=0;

	for(int i=1;i<=cat.Z;i++) if(cat.t[i]>=tt0 && cat.t[i]<=tt1) N+=1;

	integral=0.0;
	for(int i=1;i<=cat.Z;i++) dumrate[i]=rate[i]=0.0;

	err=0;

	#ifdef _CRS_MPI
		MPI_File fhw_foret1, fhw_foret2, fhw_forex, fhw_cmb;
		MPI_Status status;

		all_snapshots=0;	//fixme implement this.

		if (printall_cmb) {
			sprintf(fname, "%s.dat",printall_cmb);
			MPI_File_open(MPI_COMM_WORLD, fname,
						  MPI_MODE_CREATE|MPI_MODE_WRONLY,
						  MPI_INFO_NULL, &fhw_cmb);
		}

		if (printall_forex) {
			sprintf(fname, "%s.dat",printall_forex);
			MPI_File_open(MPI_COMM_WORLD, fname,
						  MPI_MODE_CREATE|MPI_MODE_WRONLY,
						  MPI_INFO_NULL, &fhw_forex);
		}

		if (printall_foret) {
			sprintf(fname, "%s_cumu.dat",printall_foret);
			MPI_File_open(MPI_COMM_WORLD, fname,
						  MPI_MODE_CREATE|MPI_MODE_WRONLY,
						  MPI_INFO_NULL, &fhw_foret1);

			sprintf(fname, "%s_rates.dat",printall_foret);
			MPI_File_open(MPI_COMM_WORLD, fname,
						  MPI_MODE_CREATE|MPI_MODE_WRONLY,
						  MPI_INFO_NULL, &fhw_foret2);
		}
	#else
		if(procId == 0) {
			if (printall_cmb) {
				sprintf(fname, "%s.dat",printall_cmb);
				fcmb=fopen(fname,"w");
			}
			if (printall_forex) {
				sprintf(fname, "%s.dat",printall_forex);
				fforex=fopen(fname,"w");
			}
			if (printall_foret) {
				sprintf(fname, "%s_cumu.dat",printall_foret);
				fforet1=fopen(fname,"w");
				sprintf(fname, "%s_rates.dat",printall_foret);
				fforet2=fopen(fname,"w");
			}
		}
	#endif

	if(procId == 0) {

		if (all_snapshots){
			sprintf(fname, "%s_allsnaps.dat",print_forex0);
			foutallsnap=fopen(fname,"w");
		}

		if (print_LL) {
			sprintf(fname, "%s.dat",print_LL);
			fLLev=fopen(fname,"w");
		}
		if (print_cmb) {
			sprintf(print_cmb,"%s.dat", print_cmb0);
			if (flags.afterslip) sprintf(print_cmbpost,"%s_post.dat", print_cmb0);

		}
		if (print_forex0) {
			sprintf(print_forex,"%s.dat", print_forex0);
		}
		if (print_foret) {
			sprintf(fname, "%s.dat",print_foret);
			fforet_avg=fopen(fname,"w");
		}
	}

	#ifdef _CRS_MPI
		int rootPartitionSize = 0;

		partitionSize = roundUpFrac((double)Nsur / (double)numProcs);

		if((partitionSize * numProcs) > Nsur) {
			--partitionSize;

			rootPartitionSize = partitionSize + (Nsur - (partitionSize * numProcs));

			if(procId == 0) {
				partitionSize = rootPartitionSize;

				start = 1;
			}
			else {
				start = rootPartitionSize + ((procId-1) * partitionSize) + 1;
			}
		}
		else {
			start = (procId * partitionSize) + 1;
		}

		end = start + partitionSize;

		*seed = (*seed) * (procId+numProcs);
	#else
		start = 1;
		end = Nsur + 1;
	#endif

//	// [Fahad]: FIXME -- We need to come up with a good seed generation algorithm for the MPI code
//	// [Fahad]: FIXME - For Test A2, first call ...
	long seeds[] = {-956111019, -1383064173, -25303387, -1130426989, -1321121682, -137071578, -1882507103, -1846814569, -78114812};

	for(int nsur = start; nsur < MIN(end, Nsur+1); nsur++) {
		#ifdef _CRS_MPI
//			*seed = newSeed * (long)nsur;
			// [Fahad]: FIXME -- For testing only ...
			*seed = seeds[nsur-1];
		#endif

		for (int n=1; n<=NgridT; n++) ev_x[n]=0.0;
		for(int i=1;i<=cat.Z;i++) dumrate[i]=0.0;

		if (all_snapshots){
			for (int t=0; t<Ntts; t++){
				for (int n=1; n<=NgridT; n++){
					nev_allsnapdum[t][n]=0.0;
				}
			}
		}


		print_screen("%d...",nsur);

		if (all_gammas0) gammas0= (multiple_input_gammas)? all_gammas0[nsur] : *all_gammas0;
		//if (flags.sample_all), each iteration corresponds to a focal mechanism. Otherwise, which_recfault=0 means: choose random one.
		which_recfault= flags.sample_all? nsur : 0;

		//Set starting rates:
		if(fromstart) {
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, flags,
								   times, Nm, Na, crst, AllCoeff, NTScont, focmec,
								   fmzonelim, NFM, seed, tstart, tt1, refresh && nsur==start,
								   which_recfault);

			for(int n=1; n<=NgridT; n++) {
				gammas[n]= (uniform_bg_rate)? ta/Asig : gammas0[n];
			}

			err = rate_state_evolution(cat, times, DCFSrand, DCFS, tstart, tt0, tt0-tstart, Asig, ta,
									  (int *) 0, ev_x_dum, (double *) 0, (double *) 0, (double **) 0,
									  NgridT, NTScont, Nm, gammas, crst.rate0,
									  dumrate, 1);

			for(int i=1; i<=NgridT; i++) {
				ev_x_pre[i]+=ev_x_dum[i]/(1.0*Nsur);
			}
		}
		else {
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, flags,
								   times, Nm, Na, crst, AllCoeff, NTScont, focmec,
								   fmzonelim, NFM, seed, tt0, tt1, refresh && nsur==start,
								   which_recfault);

			for(int n=1; n<=NgridT; n++) {
				gammas[n]= (uniform_bg_rate)? ta/Asig : gammas0[n];
			}
		}

		current_main=0;
		tnow=tt0;

		err+=rate_state_evolution(cat, times, DCFSrand, DCFS, tt0, tt1, dtstep,
								 Asig, ta, 0, ev_x, nev+1, rev+1, nev_allsnapdum, NgridT, NTScont,
								 Nm, gammas, crst.rate0, dumrate, 1);

		for(int t=1; t<=Ntts; t++) {
			nev_avg[t]+=nev[t]/(1.0*Nsur);
			rev_avg[t]+=rev[t]/(1.0*Nsur);

			if (all_snapshots){
				for (int m=1; m<=NgridT; m++) {
					nev_allsnapshots[t][m]+=nev_allsnapdum[t-1][m]/(1.0*Nsur);
				}
			}
		}
//		for(int i=1; i<=NgridT; i++) {
//			ev_x[i]=ev_x_dum[i];
//		}
		for(int i=1;i<=cat.Z;i++) {
			if(cat.t[i]>=tt0 && cat.t[i]<tt1) {
				rate[i]+=dumrate[i]/(1.0*Nsur);
			}
		}

		if(err==1) break;

		for (int i=1; i<=NgridT; i++) {
			ev_x_avg[i]+=ev_x[i]/(1.0*Nsur);
		}

		if(print_cmb || printall_cmb) {
			sum_DCFS(DCFS, &cmb, Nm, NgridT);
			if (flags.afterslip) sum_DCFSrand(DCFSrand, &cmbpost, NTScont, NgridT);
			if(print_cmb) {
				for (int i=1; i<=NgridT; i++) {
					cmb_avg[i]+=cmb[i]*(1.0/Nsur);
					if (flags.afterslip) cmbpost_avg[i]+=cmbpost[i]*(1.0/Nsur);
				}
			}

			if(printall_cmb) {
				convert_geometry(crst,cmb, &ev_x_new, 0, 0);	//convert to output geometry
				#ifdef _CRS_MPI
					int offset = ((nsur-1)*(NgridT_out)*sizeof(double));

					// [Fahad]: Buffer indices have a '+1' because these are nrutils dvector types.
					MPI_File_write_at(fhw_cmb, offset, ev_x_new+1, NgridT_out,
									  MPI_DOUBLE, &status);
				#else
					for (int n=1; n<=NgridT_out; n++) fprintf(fcmb, "%.6e\t", ev_x_new[n]);
					if (nsur <Nsur) fprintf(fcmb, "\n");
					fflush(fcmb);
				#endif
			}
		}

		if (printall_forex) {
			convert_geometry(crst, ev_x, &ev_x_new, 1, 0);	//convert to output geometry
			#ifdef _CRS_MPI
				for(int n = 1; n <= NgridT_out; n++) {
					ev_x_new[n] *= r0/NgridT;
				}

				int offset = ((nsur-1)*(NgridT_out)*sizeof(double));

				// [Fahad]: Buffer indices have a '+1' because these are nrutils dvector types.
				MPI_File_write_at(fhw_forex, offset, ev_x_new+1, NgridT_out,
								  MPI_DOUBLE, &status);
			#else
				for (int n=1; n<=NgridT_out; n++) fprintf(fforex, "%.6e\t", ev_x_new[n]*r0/NgridT);
				if (nsur <Nsur) fprintf(fforex, "\n");
				fflush(fforex);
			#endif
		}

		if (printall_foret) {
			#ifdef _CRS_MPI
				for (int t=1; t<=Ntts; t++) {
					nev[t] *= r0 / NgridT;
					rev[t] *= r0 / NgridT;
				}

				int offset = ((nsur-1)*(Ntts)*sizeof(double));

				// [Fahad]: Buffer indices have a '+1' because these are nrutils dvector types.
				MPI_File_write_at(fhw_foret1, offset, nev+1, Ntts,
								  MPI_DOUBLE, &status);
				MPI_File_write_at(fhw_foret2, offset, rev+1, Ntts,
								  MPI_DOUBLE, &status);
			#else
				for (int t=1; t<=Ntts; t++) {
					fprintf(fforet1, "%lf\t",nev[t]*r0/NgridT);
					fprintf(fforet2, "%lf\t",rev[t]*r0/NgridT);
				}
				if (nsur <Nsur) fprintf(fforet1, "\n");
				if (nsur <Nsur) fprintf(fforet2, "\n");
				fflush(fforet1);
				fflush(fforet2);
			#endif
		}
	}

	#ifdef _CRS_MPI
		double *recv_rate, *recv_nev_avg, *recv_rev_avg, *recv_cmb_avg, *recv_ev_x_avg;
		recv_rate = dvector(1, cat.Z);
		recv_nev_avg = dvector(1, Ntts);
		recv_rev_avg = dvector(1, Ntts);
		recv_cmb_avg = dvector(1,NgridT);
		recv_ev_x_avg = dvector(1, NgridT);

		MPI_Allreduce(rate, recv_rate, cat.Z+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(nev_avg, recv_nev_avg, Ntts+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(rev_avg, recv_rev_avg, Ntts+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(cmb_avg, recv_cmb_avg, NgridT+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(ev_x_avg, recv_ev_x_avg, NgridT+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		free_dvector(rate, 1, cat.Z);
		rate = recv_rate;
		free_dvector(nev_avg, 1, Ntts);
		nev_avg = recv_nev_avg;
		free_dvector(rev_avg, 1, Ntts);
		rev_avg = recv_rev_avg;
		free_dvector(cmb_avg, 1, NgridT);
		cmb_avg = recv_cmb_avg;
		free_dvector(ev_x_avg, 1, NgridT);
		ev_x_avg = recv_ev_x_avg;
	#endif

	//calculate average rate and LL:
	if(!err) {
		if (printall_cmb) {
			#ifdef _CRS_MPI
				MPI_File_close(&fhw_forex);
			#else
				fclose(fcmb);
			#endif
		}

		if (printall_forex) {
			#ifdef _CRS_MPI
				MPI_File_close(&fhw_forex);
			#else
				fclose(fforex);
			#endif
		}

		if (printall_foret) {
			#ifdef _CRS_MPI
				MPI_File_close(&fhw_foret1);
				MPI_File_close(&fhw_foret2);
			#else
				fclose(fforet1);
				fclose(fforet2);
			#endif
		}

		if(procId == 0) {
			if (print_foret) {
				for (int t=1; t<=Ntts; t++) fprintf(fforet_avg, "%lf\t%lf\t%lf\n",tt0+dtstep*t,nev_avg[t]*r0/(1.0*NgridT),rev_avg[t]*r0/(1.0*NgridT));
				fclose(fforet_avg);
			}
		}

		if (all_snapshots) {

			for(int t=1; t<=Ntts; t++) {
				convert_geometry(crst, nev_allsnapshots[t], &ev_x_new, 1, 0);
				if (flatten_allsnaps) {
					flatten_outgrid(crst, ev_x_new, &flat_grid, &NgridTsnaps);
				}
				else {
					flat_grid=ev_x_new;
					NgridTsnaps= NgridT_out;
				}
				if(procId == 0) {
					for(int n=1; n<=NgridTsnaps; n++) {
						fprintf(foutallsnap, "%.5e\t", flat_grid[n]*=r0/NgridT);
					}
					fprintf(foutallsnap, "\n");
				}
			}
		}


		if (print_forex) {
			convert_geometry(crst, ev_x_avg, &ev_x_new, 1, 0);
			for(int n=1; n<=NgridT_out; n++) {
				ev_x_new[n]*=r0/NgridT;
			}
			if(procId == 0) {
				csep_forecast(print_forex, crst, ev_x_new, 0);
			}
		}

		if (print_cmb) {
			convert_geometry(crst, cmb_avg, &ev_x_new, 0, 0);
			if(procId == 0) {
				csep_cmbmap(print_cmb, crst, ev_x_new, 0);
				if (flags.afterslip) csep_cmbmap(print_cmbpost, crst, cmbpost_avg, 0);
			}
		}
		if (print_LL || LL){
			Ldum=0.0;

			for(int j=1;j<=cat.Z;j++) if(cat.t[j]>=tt0 && cat.t[j]<tt1) {
				Ldum+=log(r0*rate[j]);
				if(procId == 0) {
					if (print_LL) fprintf(fLLev,"%.10e \t%.5e \t%.5e \t%.5e \t %.5e\n",cat.t[j], cat.lat0[j], cat.lon0[j], cat.depths0[j], log(r0*rate[j]));
				}
			}

			if(procId == 0) {
				if (print_LL) fclose(fLLev);
			}

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
	free_dvector(ev_x,1,NgridT);
	free_dvector(ev_x_avg,1,NgridT);
	free_dvector(ev_x_new,1,NgridT_out);
	free_dvector(cmb,1,NgridT);
	free_dvector(cmbpost,1,NgridT);
	free_dvector(cmb_avg,1,NgridT);
	free_dvector(nev,1,Ntts);
	free_dvector(rev,1,Ntts);
	free_dvector(nev_avg,1,Ntts);
	free_dvector(rev_avg,1,Ntts);

	return(err);
}

int CRSLogLikelihood(double *LL, double *Ldum0_out, double *Nev, double *I, double *r_out,
					 int Nsur, struct pscmp *DCFS, struct eqkfm *eqkfm_aft,
					 struct eqkfm *eqkfm0, struct flags flags,
					 struct crust crst, struct Coeff_LinkList *AllCoeff, int NTScont,
					 int Nm, int Na, int NgridT, double **focmec, int *fmzonelim, int NFM,
					 long *seed, struct catalog cat, double *times, double tstart, double tt0,
					 double tt1, double tw, double Mag_main, double Asig, double ta, double r0, int fixr,
					 double *gammas0, double **all_new_gammas, int fromstart, char * printall_cmb,
					 char *printall_forex, int refresh) {

	//recfault= [0,1,2] means: don't vary rec. fault, vary (choose random one), vary and sample all catalog in order.
	//tstart=starting time of entire period considered (may be <tt0). only used if fromstart==1;
	//if multiple_input_gammas==1, expect allgammas0[1...Nsur][1...NgridT]; else, allgammas0is pointer to vector: [1...NgridT].
	//same for multiple_input_gammas (if multiple_output_gammas, return gammas that would give average rate).
	//gammas0 can be an array of starting values for gamma; they will be non steady-state values. Otherwise, ...

	/*INPUT:-------------------------------------------------------------//
	 *
	 * Nsur= no. of iterations
	 * DCFS= array containing stress steps.
	 * eqkfm_aft= array of afterslip:
	 * 		if splines==1
	 * 		if splines==0
	 * eqkfm0= array of mainshocks;
	 * crst= general model info.
	 * AllCoeff= Okada Coefficients for all mainshock slip models;
	 * NTScont=no. of continuous time steps (size of times2)
	 * Nm=no. of mainshocks;	size of eqkfm0
	 * Na= no. of events with afterslip
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

	// [Fahad] Variables used for MPI.
	int procId = 0, numProcs = 1;
	int start, end, partitionSize;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	#endif

	static double **DCFSrand;
	static double *dumrate, *gammas, *rate, *ev_x;
	double sum, sum1, sum2, integral, Ldum0, r;
	static int first_timein=1, N;
	double Ntot, Itot, LLdum0tot;
	int err;
	int current_main, j0, which_recfault;
	double tnow;
	FILE *fforex, *fcmb;

	#ifdef _CRS_MPI
		// [Fahad]: The dmatrix 'all_new_gammas' has to be linearized for the
		//		  : MPI communication routine that consolidates values from
		//		  : all ranks.
		size_t allNewGammasSize;
		size_t localNewGammasSize;

		double* linearizedAllNewGammas;
		double* linearizedLocalNewGammas;
	#endif

	if (LL) print_screen("Calculating LL for Asig=%lf, ta=%lf ...", Asig, ta);

	if (first_timein==1){
		print_logfile("Setting up variables in CRSLogLikelihood...\n");
		#ifndef _CRS_MPI
			// [Fahad] Assignment happens at the end of the function when MPI is being used.
			first_timein=0;
		#endif
		DCFSrand= (flags.afterslip) ? dmatrix(0,NTScont,1,NgridT) : NULL;
		dumrate=dvector(1,cat.Z);
		rate=dvector(1,cat.Z);
		gammas=dvector(1,NgridT);
		ev_x=dvector(1,NgridT);
	}

	sum=sum1=sum2=0.0;
	integral=0.0;
	for(int i=1;i<=cat.Z;i++) {
		rate[i]=0.0;
	}

	sum=sum1=sum2=0;
	err=0;

	if (printall_cmb) fcmb=fopen(printall_cmb,"w");
	if (printall_forex) fforex=fopen(printall_forex,"w");

	#ifdef _CRS_MPI
		int rootPartitionSize = 0;

		if(first_timein != 1) {
			partitionSize = roundUpFrac((double)Nsur / (double)numProcs);

			if((partitionSize * numProcs) > Nsur) {
				--partitionSize;

				rootPartitionSize = partitionSize + (Nsur - (partitionSize * numProcs));

				if(procId == 0) {
					partitionSize = rootPartitionSize;

					start = 1;
				}
				else {
					start = rootPartitionSize + ((procId-1) * partitionSize) + 1;
				}
			}
			else {
				start = (procId * partitionSize) + 1;
			}

			end = start + partitionSize;
		}
		else {
			start = 1;
			end = Nsur + 1;
		}

		if(all_new_gammas) {
			if(rootPartitionSize) {
				allNewGammasSize   = (rootPartitionSize * numProcs) * NgridT;
				localNewGammasSize = rootPartitionSize * NgridT;
			}
			else {
				allNewGammasSize   = Nsur * NgridT;
				localNewGammasSize = (MIN(end, Nsur+1) - start) * NgridT;
			}

			linearizedAllNewGammas   = (double*) malloc(allNewGammasSize   * sizeof(double));
			linearizedLocalNewGammas = (double*) malloc(localNewGammasSize * sizeof(double));
		}

		*seed = (*seed) * (procId+numProcs);
	#else
		start = 1;
		end = Nsur + 1;
	#endif

//	// [Fahad]: FIXME -- We need to come up with a good seed generation algorithm for the MPI code
//	// [Fahad]: FIXME - For Test A2, first call ...
	long seeds[] = {-956111019, -1383064173, -25303387, -1130426989, -1321121682, -137071578, -1882507103, -1846814569, -78114812};

	for(int nsur = start; nsur < MIN(end, Nsur+1); nsur++) {
		#ifdef _CRS_MPI
			// [Fahad]: FIXME -- For testing only ...
			if(first_timein != 1) {
				*seed = seeds[nsur-1];
			}
		#endif


		//if (flags.sample_all), each iteration corresponds to a focal mechanism. Otherwise, which_recfault=0 means: choose random one.
		which_recfault= flags.sample_all? nsur : 0;	//which_recfault=0 means: choose random one.

		for(int i=1;i<=cat.Z;i++) dumrate[i]=0.0;

		//Set starting rates:
		if (fromstart){
		// todo [coverage] this block is never tested
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, flags,
								   times, Nm, Na, crst, AllCoeff, NTScont, focmec,
								   fmzonelim, NFM, seed, tstart, tt1,
								   refresh && nsur==start /*&& first_timein*/, which_recfault);

			for(int n=1; n<=NgridT; n++) {
				gammas[n]= (gammas0)? gammas0[n] : ta/Asig;	//if gammas0 NULL, use uniform background rate (steady state).
			}

			// the time step passed to rate_state_evolution (dt_step) must be larger than tt1-tt0 to avoid having 2 time steps;
			// tt0-start may not be enough due to floating point error --> use tt0-start+1.0. Also done below.
			err = rate_state_evolution(cat, times, DCFSrand, DCFS, tstart, tt0, tt0-start+1.0, Asig, ta,
									  (int *) 0, (double *) 0, (double *) 0, (double *) 0, (double **) 0,
									  NgridT, NTScont, Nm, gammas, (double *) 0,
									  dumrate, 1);
		}
		else{
			calculateDCFSperturbed(DCFSrand, DCFS, eqkfm_aft, eqkfm0, flags,
								   times, Nm, Na, crst, AllCoeff, NTScont, focmec,
								   fmzonelim, NFM, seed, tt0, tt1,
								   refresh && nsur==start /*&& first_timein*/, which_recfault);

			for(int n=1; n<=NgridT; n++) {
				gammas[n]= (gammas0)? gammas0[n] : ta/Asig;	//if gammas0 NULL, use uniform background rate (steady state).
			}
		}

		//Calculate seismicity evolution (skipping a time window after each mainshock):
		current_main=0;
		tnow=tt0;
		//find first event which causes a time window of incomplete seismicity:
		while (current_main<Nm && eqkfm0[current_main].t<tt0 && eqkfm0[current_main].mag<Mag_main) current_main++;
		while (current_main<Nm && eqkfm0[current_main].t<tt1){
			if (tnow<eqkfm0[current_main].t){
				//evolve seismicity up to next large event:
				err += rate_state_evolution(cat, times, DCFSrand, DCFS, tnow, eqkfm0[current_main].t, eqkfm0[current_main].t-tnow+1.0,
										   Asig, ta, 0, 0, &sum, 0, NULL, NgridT, NTScont, Nm, gammas,
										   crst.rate0, dumrate, 1);
				integral += (sum)/(1.0*Nsur);

				//evolve seismicity during a time window tw:
				err += rate_state_evolution(cat, times, DCFSrand, DCFS, eqkfm0[current_main].t, eqkfm0[current_main].t+tw,
										tw+1.0, Asig, ta, 0, 0, &sum, 0, NULL,
									   NgridT, NTScont, Nm, gammas, crst.rate0,
									   dumrate, 1);
				tnow=eqkfm0[current_main].t+tw;
			}
			else {
				//the condition below will be true if the current large event is still within the tw of the previous one.
				//(Actually, if tw has a fixed value this is always the case).
				if (tnow<eqkfm0[current_main].t+tw){
					err += rate_state_evolution(cat, times, DCFSrand, DCFS, tnow, eqkfm0[current_main].t+tw, eqkfm0[current_main].t+tw-tnow+1.0,
							Asig, ta, 0, ev_x, &sum, 0, NULL, NgridT, NTScont, Nm, gammas, crst.rate0,
							dumrate, 1);
					tnow=eqkfm0[current_main].t+tw;
				}
			}

			//find next event large enough to causes a time window of incomplete seismicity:
			current_main+=1;
			while (current_main<Nm && eqkfm0[current_main].mag<Mag_main) current_main++;
			if (err) break;
		}
		if (tnow<tt1){
			err += rate_state_evolution(cat, times, DCFSrand, DCFS, tnow, tt1, tt1-tnow+1.0, Asig, ta, 0,
									   ev_x, &sum, 0, NULL, NgridT, NTScont, Nm, gammas,
									   crst.rate0, dumrate, 1);

			integral+=(sum)/(1.0*Nsur);
		}

		if (err==1) break;

		for(int i=1;i<=cat.Z;i++) {
			if(cat.t[i]>=tt0 && cat.t[i]<tt1) {
				rate[i]+=1.0*dumrate[i]/(1.0*Nsur);
			}
		}

		if (all_new_gammas) {
			for (int n=1; n<=NgridT; n++) {
				#ifdef _CRS_MPI
					linearizedLocalNewGammas[(nsur-start)*NgridT + (n-1)] = gammas[n];
				#else
					all_new_gammas[nsur][n]=gammas[n];
				#endif
			}
		}

		// todo [coverage] this block (2x) is never tested
		if (printall_cmb) {
			for (int n=1; n<=NgridT; n++) fprintf(fcmb, "%lf\t", DCFS[0].cmb[n]);
			if (nsur <Nsur) fprintf(fcmb, "\n");
		}
		if (printall_forex) {
			for (int n=1; n<=NgridT; n++) fprintf(fforex, "%lf\t", ev_x[n]);
			if (nsur <Nsur) fprintf(fforex, "\n");
		}
	}

	#ifdef _CRS_MPI
		double temp_integral;
		double *recv_rate, *recv_rates_x;

		recv_rate = dvector(1, cat.Z);

		MPI_Allreduce(&integral, &temp_integral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(rate, recv_rate, cat.Z+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		integral = temp_integral;
		free_dvector(rate, 1, cat.Z);
		rate = recv_rate;

		// [Fahad]: TODO -- This is way too complicated. Try to come up with a more
		//				 -- elegant algorithm. Add comments in the meanwhile ...
		if (all_new_gammas) {
			MPI_Allgather(linearizedLocalNewGammas, localNewGammasSize, MPI_DOUBLE,
						  linearizedAllNewGammas,   localNewGammasSize, MPI_DOUBLE,
						  MPI_COMM_WORLD);

			// Copy values from the linearized matrix to dmatrix
			if(rootPartitionSize) {	// If root has a larger partition size than the other ranks
				// We need to skip empty portions of NgridT size in the linearized array,
				// for all non-root ranks.
				for(int i = 0; i < numProcs; ++i) {
					int numRows, rankIndexLinear, rankIndexMatrix;

					if(i == 0) {
						numRows = rootPartitionSize;
					}
					else {
						numRows = roundUpFrac((double)Nsur / (double)numProcs) - 1;
					}

					rankIndexLinear = i * rootPartitionSize * NgridT;
					rankIndexMatrix = rootPartitionSize + ((i-1) * numRows);

					for(int j = 0; j < numRows; ++j) {
						for(int k = 0; k < NgridT; ++k) {
							all_new_gammas[rankIndexMatrix + j+1][k+1] = linearizedAllNewGammas[rankIndexLinear + (j*NgridT) + k];
						}
					}
				}
			}
			else {
				for(int i = 0; i < Nsur; ++i) {
					for(int j = 0; j < NgridT; ++j) {
						all_new_gammas[i+1][j+1] = linearizedAllNewGammas[(i*NgridT) + j];
					}

				}
			}

			free(linearizedLocalNewGammas);
			free(linearizedAllNewGammas);
		}
	#endif

	//calculate average rate and LL:
	if (!err){
		if (printall_cmb) fclose(fcmb);
		if (printall_forex) fclose(fforex);

		Ldum0=0.0;
		N=0;
		current_main=0;
		tnow=tt0;
		j0=1;

		FILE *fout;
		if(procId == 0) {
			fout = fopen("LLs.dat","w");
		}

		while (current_main<Nm && eqkfm0[current_main].t<tt0 && eqkfm0[current_main].mag<Mag_main) current_main++;
		while (current_main<Nm && eqkfm0[current_main].t<tt1){
			if (tnow<=eqkfm0[current_main].t){

				while(j0<=cat.Z && cat.t[j0]<=eqkfm0[current_main].t){
					if(cat.t[j0]>=tnow) {
						N+=1;
						Ldum0+=log(rate[j0]);
						if(procId == 0) {
							fprintf(fout,"%lf\t%lf\t%10e\n", cat.t[j0], cat.mag[j0], rate[j0]);
						}
					}
					j0+=1;
				}
			}
			tnow=eqkfm0[current_main].t+tw;
			current_main+=1;
			while (current_main<Nm && eqkfm0[current_main].mag<Mag_main) current_main++;
		}

		if (tnow<tt1){
			for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<tt1) {
				N+=1;
				Ldum0+=log(rate[j]);
				if(procId == 0) {
					fprintf(fout,"%lf\t%lf\t%10e\n", cat.t[j],cat.mag[j], rate[j]);
				}
			}
		}

		if(procId == 0) {
			fclose(fout);
		}

		Ntot= (Nev) ? N+*Nev : N;
		Itot= (I)? integral/NgridT + *I : integral/NgridT;
		LLdum0tot= (Ldum0_out) ? Ldum0+*Ldum0_out : Ldum0;
		r = (fixr)? r0 : Ntot/Itot;	// Ntot/Itot can be derived analytically (gives maximum of LL wrt r).

		if(r==0.0) {
			print_screen("ERROR: ta=%lf  Asig=%lf r=%e\n",ta,Asig,r);
			print_logfile("ERROR: ta=%lf  Asig=%lf r=%e\n",ta,Asig,r);
		}

		if (LL) *LL=LLdum0tot+ Ntot*log(r) - r*Itot;
		if (Ldum0_out) *Ldum0_out=LLdum0tot;
		if (I) *I=Itot;
		if (Nev) *Nev=Ntot;
		if (r_out) *r_out=r;

		if (LL) print_screen("LL=%lf", *LL);
	}

	print_screen("\n");
	if (first_timein) print_logfile("done.\n");

	#ifdef _CRS_MPI
		first_timein = 0;
	#endif

	return(err);
}
