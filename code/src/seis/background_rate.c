/*
 * background_rate.c
 *
 *  Created on: Nov 12, 2013
 *      Author: camcat
 */

#include "background_rate.h"

int background_rate(char *catfile, struct crust *crst_in, struct tm reftime, double Mcut, double Mmain,  double t0, double t1, double dR, double dZ, int ord){
	/* by convention, if Mcut>=20 the cutoff magnitude will be calculated based on completeness.
	 * ord=1,2 indicatest if closest or second closest event should be used for smoothing distance.
	 * does not decluster catalog.
	 */

	struct catalog cat;
	struct crust crst=*crst_in;
	double *zlist;
	int NP=crst.nLat*crst.nLon;
	double T, *rate_h, *rate_v;
	double Mc, b, Mc_final=crst.mags[1]+0.5*(crst.dmags);
	int ne, h_ind, v_ind;
	cat.Mc=Mcut;

	zlist=dvector(1,crst.nD);
	for (int i=1; i<=crst.nD; i++) zlist[i]=crst.depth[1+(i-1)*NP];

	readZMAP(&cat, NULL, NULL, catfile, crst, reftime, t0, t1, t0, t1, 10, 0.0, dR, dZ, 0.0, 0);
	if (cat.Z==0) return 1;

	//avoid events following the first large earthquake.
	for (int e=1; e<=cat.Z; e++) {
		if (cat.mag[e]>=Mmain){
			cat.tend=cat.t[e];
			cat.Z=e;
		}
	}
	T=cat.tend-cat.tstart;

	if (cat.Z<=100){
		if (cat.Z<=1 || T<1.0)	{
			if (verbose_level) printf("**Warning: to few events, or too short time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
			if (flog) {
				fprintf(flog,"**Warning: to few events, or too short a valid time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
				fflush(flog);
			}
			return 1;
		}
		else{
			if (verbose_level || flog){
				if (verbose_level) printf("**Warning: only %ld events, used to calculate background seismicity rate.** \n", cat.Z);
				if (flog) {
					fprintf(flog,"**Warning: only %ld events, used to calculate background seismicity rate.** \n", cat.Z);
					fflush(flog);
				}
			}
		}
	}

	//find rate by counting events of required magnitude; is too few (<10), use GR instead.

	(*crst_in).r0=0.0;
	for (int e=1; e<=cat.Z; e++) if (cat.mag[e]>=Mc_final) (*crst_in).r0+=1;
	if ((*crst_in).r0>=10) (*crst_in).r0*=1.0/T;
	else{
		if (Mcut>=20) {
			(*crst_in).r0=((double) cat.Z-1.0)/T;
			b=cat.b;
			Mc=cat.Mc;
		}
		else {
			Mc=Mc_maxcurv(cat.mag+1, cat.Z)+0.2;
			b=calculatebvalue(cat.mag+1, cat.Z, Mc);
			ne=0;
			for (int e=1; e<=cat.Z; e++) if (cat.mag[e]>=Mc) ne++;
			(*crst_in).r0=((double) ne)/T;
		}
		(*crst_in).r0*=pow(10,b*(Mc-Mc_final));
	}


	rate_h=Helmstetter_cat(cat, crst, NULL, ord);
	rate_v=fit_depth(zlist, zlist[2]-zlist[1], crst.nD, cat.depths0, cat.verr, NULL, cat.Z);

	normv(rate_h, NP);
	normv(rate_v, crst.nD);

	(*crst_in).rate0=dvector(1,crst.N_allP);
	for (int p=1; p<=crst.N_allP; p++) {
		h_ind=(p-1)%NP+1;
		v_ind=(p-1)/NP+1;
		(*crst_in).rate0[p]=rate_h[h_ind]*rate_v[v_ind]*crst.N_allP;
	}

	free_cat(cat);
	free_dvector(zlist,1,crst.nD);
	free_dvector(rate_h, 1, NP);
	free_dvector(rate_v, 1, crst.nD);

	return 0;

}

int background_rate2(char *catfile, struct crust *crst_in, struct tm reftime, double Mcut, double Mmain, double *target_mags, double *target_rates, int Ntarget, double t0, double t1, double dR, double dZ, double min_smoothing, int ord){
	/* by convention, if Mcut>=20 the cutoff magnitude will be calculated based on completeness.
	 * target_mags= array containing magnitudes above which rates should be estimated (counting events above Mw, or using GR if there are too few). [0....Ntarget-1];
	 * target_rates will contain the results. If NULL (or if target_mags=NULL), ignored.
	 * ord=1,2 indicatest if closest or second closest event should be used for smoothing distance.
	 * declusters catalog using KG window method, and scales remaining events to compensate for removing background seismicity.
	 */

	struct catalog cat;
	struct crust crst=*crst_in;
	double *zlist;
	int NP=crst.nLat*crst.nLon;
	double T, *rate_h, *rate_v;
	double Mc, b, Mc_final=crst.mags[1]-0.5*(crst.dmags), rcat;
	int ne, h_ind, v_ind, *sel, sel_no=0;
	double *weights=NULL;
	cat.Mc=Mcut;

	zlist=dvector(1,crst.nD);
	for (int i=1; i<=crst.nD; i++) zlist[i]=crst.depth[1+(i-1)*NP];

	//todo remove first line.
	if (countcol(catfile)==8) read_RS(catfile, &cat, crst, -2, t0, 0.0, 0.0, t1, 0.0, (struct eqkfm **) 0, 0, NULL, 0);
	else readZMAP(&cat, (struct eqkfm **) 0, NULL, catfile, crst, reftime, t0, t1, t0, t1, 10, 0.0, dR, dZ, 0.0, 0);
	if (cat.Z==0) return 1;

	sel=decluster_catalog(cat, Mmain, &weights, 0);

	for (int i=1; i<=cat.Z; i++) {
		cat.err[i]=fmax(cat.err[i], min_smoothing);	//values used later (in Helmstetter.c).
		sel_no+=sel[i];
	}

	T=cat.tend-cat.tstart;
	if (sel_no<=100){
		if (sel_no<=1 || T<1.0)	{
			if (verbose_level) printf("**Warning: to few events, or too short time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
			if (flog) {
				fprintf(flog,"**Warning: to few events, or too short a valid time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
				fflush(flog);
			}
			return 1;
		}
		else{
			if (verbose_level || flog){
				if (verbose_level) printf("**Warning: only %ld events used to calculate background seismicity rate.** \n", cat.Z);
				if (flog) {
					fprintf(flog,"**Warning: only %ld events used to calculate background seismicity rate.** \n", cat.Z);
					fflush(flog);
				}
			}
		}
	}

	//find rate of events in catalog:
	if (Mcut>=20) {
		rcat=((double) cat.Z-1.0)/T;
		b=cat.b;
		Mc=cat.Mc;
	}
	else {
		Mc=Mc_maxcurv(cat.mag+1, cat.Z)+0.2;
		b=calculatebvalue(cat.mag+1, cat.Z, Mc);
		ne=0;
		for (int e=1; e<=cat.Z; e++) if (cat.mag[e]>=Mc) ne++;
		rcat=((double) ne)/T;
	}
	rcat*=pow(10,b*Mc);

	//find rate by counting events of required magnitude; is too few (<10), use GR instead.
	if (target_rates && target_mags){

		for (int i=0; i<Ntarget; i++){
			target_rates[i]=0.0;
			for (int e=1; e<=cat.Z; e++) if (cat.mag[e]>=target_mags[i]) target_rates[i]+=1;
			if (target_rates[i]>=10) target_rates[i]*=1.0/T;
			else (*crst_in).r0*=rcat*pow(10,-b*target_mags[i]);
		}
	}

	//do same for crst:
	//too few events (<10), use GR instead.
	(*crst_in).r0=0.0;
	for (int e=1; e<=cat.Z; e++) if (cat.mag[e]>=Mc_final) (*crst_in).r0+=1;
	if ((*crst_in).r0>=10) (*crst_in).r0*=1.0/T;
	else (*crst_in).r0*=pow(10,-b*Mc_final);

	rate_h=Helmstetter_cat(cat, crst, weights, ord);
	rate_v=fit_depth(zlist, zlist[2]-zlist[1], crst.nD, cat.depths0, cat.verr, weights, cat.Z);

	normv(rate_h, NP);
	normv(rate_v, crst.nD);

	(*crst_in).rate0=dvector(1,crst.N_allP);
	for (int p=1; p<=crst.N_allP; p++) {
		h_ind=(p-1)%NP+1;
		v_ind=(p-1)/NP+1;
		(*crst_in).rate0[p]=rate_h[h_ind]*rate_v[v_ind]*crst.N_allP;
	}

	free_cat(cat);
	free_dvector(zlist,1,crst.nD);
	free_dvector(rate_h, 1, NP);
	free_dvector(rate_v, 1, crst.nD);

	return 0;

}

int background_rate3(char *catfile, struct crust *crst_in, struct tm reftime, double Mcut, double Mmain,  double t0, double t1, double dR, double dZ, int ord){
	/* by convention, if Mcut>=20 the cutoff magnitude will be calculated based on completeness.
	 * ord=1,2 indicatest if closest or second closest event should be used for smoothing distance.
	 * does not decluster catalog, but only uses time period following last large (M>Mmain) event (skipping a time window given by WG algorithm).
	 * Mcut used to calculate smoothed seismicity, Mcutfinal (from crst_in) to calculate no. of events (instead of using Gutenberg Richter as in other functions.
	 */

	struct catalog cat;
	struct crust crst=*crst_in;
	double *zlist;
	int NP=crst.nLat*crst.nLon;
	double T, *rate_h, *rate_v, Tw;
	double b, Mc, ne;
	double Mc_final=crst.mags[1]+0.5*(crst.dmags);
	int h_ind, v_ind, last_main, first_out;
	double *weights=NULL;
	cat.Mc=Mcut;

	zlist=dvector(1,crst.nD);
	for (int i=1; i<=crst.nD; i++) zlist[i]=crst.depth[1+(i-1)*NP];

	readZMAP(&cat, NULL, NULL, catfile, crst, reftime, t0, t1, t0, t1, 10, 0.0, dR, dZ, 0.0, 0);
	if (cat.Z==0) return 1;

	last_main=cat.Z-1;
	while (last_main>0 && cat.t[last_main]>=t1) last_main--;
	while (last_main>0 && cat.t[last_main]>=t0 && cat.mag[last_main]<Mmain) last_main--;

	first_out=last_main;
	if (last_main>0) {
		KG74(cat.mag[last_main], NULL, &Tw);
		cat.tstart=cat.t[last_main]+Tw;
		while (first_out<cat.Z && cat.t[first_out]<cat.tstart) first_out++;
		shift_cat(&cat, first_out);
	}

	T=cat.tend-cat.tstart;
	if (cat.Z<=100){
		if (cat.Z<=1 || T<1.0)	{
			if (verbose_level) printf("**Warning: to few events, or too short time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
			if (flog) {
				fprintf(flog,"**Warning: to few events, or too short a valid time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
				fflush(flog);
			}
			return 1;
		}
		else{
			if (verbose_level || flog){
				if (verbose_level) printf("**Warning: only %ld events, used to calculate background seismicity rate.** \n", cat.Z);
				if (flog) {
					fprintf(flog,"**Warning: only %ld events, used to calculate background seismicity rate.** \n", cat.Z);
					fflush(flog);
				}
			}
		}
	}

	//find rate by counting events of required magnitude; is too few (<10), use GR instead.
	(*crst_in).r0=0.0;
	for (int e=1; e<=cat.Z; e++) if (cat.mag[e]>=Mc_final) (*crst_in).r0+=1;
	if ((*crst_in).r0>=10) (*crst_in).r0*=1.0/T;
	else{
		if (Mcut>=20) {
			(*crst_in).r0=((double) cat.Z-1.0)/T;
			b=cat.b;
			Mc=cat.Mc;
		}
		else {
			Mc=Mc_maxcurv(cat.mag+1, cat.Z)+0.2;
			b=calculatebvalue(cat.mag+1, cat.Z, Mc);
			ne=0;
			for (int e=1; e<=cat.Z; e++) if (cat.mag[e]>=Mc) ne++;
			(*crst_in).r0=((double) ne)/T;
		}
		(*crst_in).r0*=pow(10,b*(Mc-Mc_final));
	}

	rate_h=Helmstetter_cat(cat, crst, weights, ord);
	rate_v=fit_depth(zlist, zlist[2]-zlist[1], crst.nD, cat.depths0, cat.verr, weights, cat.Z);

	normv(rate_h, NP);
	normv(rate_v, crst.nD);

	(*crst_in).rate0=dvector(1,crst.N_allP);
	for (int p=1; p<=crst.N_allP; p++) {
		h_ind=(p-1)%NP+1;
		v_ind=(p-1)/NP+1;
		(*crst_in).rate0[p]=rate_h[h_ind]*rate_v[v_ind]*crst.N_allP;
	}

	if (first_out) shift_cat(&cat, 2-1*first_out);
	free_cat(cat);
	free_dvector(zlist,1,crst.nD);
	free_dvector(rate_h, 1, NP);
	free_dvector(rate_v, 1, crst.nD);

	return 0;

}

