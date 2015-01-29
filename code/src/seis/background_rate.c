/*
 * background_rate.c
 *
 *  Created on: Nov 12, 2013
 *      Author: camcat
 */

#include "background_rate.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int background_rate2(char *catfile, struct crust *crst_in, struct tm reftime,
		double Mmain, double *minmag, double *rate, double **rate_grid, double dR, double dZ, double min_smoothing, int ord) {
	/*
	 *
	 *
	 *
	 *
	 *
	 * by convention, if Mcut>=20 the cutoff magnitude will be calculated based on completeness.
	 * target_mags= array containing magnitudes above which rates should be estimated (counting events above Mw, or using GR if there are too few). [0....Ntarget-1];
	 * target_rates will contain the results. If NULL (or if target_mags=NULL), ignored.
	 * ord=1,2 indicatest if closest or second closest event should be used for smoothing distance.
	 * declusters catalog using KG window method, and scales remaining events to compensate for removing background seismicity.
	 */

	int procId = 0;
	int fileError=0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct catalog cat;
	struct crust crst=*crst_in;
	double *zlist=NULL;
	double t0;
	int NP= crst.uniform? (crst.nLat*crst.nLon) : crst.N_allP;
	double T, *rate_h=NULL, *rate_v=NULL;
	int ne, h_ind, v_ind, *sel, sel_no=0, err;
	double *weights=NULL;
	cat.Mc=20.0;	//by convention, this value will make readZMAP find completeness magnitude.

	read_firstlineZMAP(catfile, reftime, &t0);	//just to know starting time
	err=readZMAP(&cat, (struct eqkfm **) 0, NULL, catfile, crst, reftime, t0, 0.0, t0, 0.0, 10, 0.0, dR, dZ, 0.0, 0);

	if(err || cat.Z==0) {
		return 1;
	}

	//decluster catalog:
	sel=decluster_catalog(cat, Mmain, &weights, 0);

	for (int i=1; i<=cat.Z; i++) {
		cat.err[i]=fmax(cat.err[i], min_smoothing);	//values used later (in Helmstetter.c).
		sel_no+=sel[i];
	}

	T=cat.tend-cat.tstart;
	if (sel_no<=100){
		if (sel_no<=1 || T<1.0)	{
			print_screen("**Warning: to few events, or too short time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
			print_logfile("**Warning: to few events, or too short a valid time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
			return 1;
		}
		else{
			print_screen("**Warning: only %ld events used to calculate background seismicity rate.** \n", cat.Z);
			print_logfile("**Warning: only %ld events used to calculate background seismicity rate.** \n", cat.Z);
		}
	}
	else {
		print_logfile("**Background seismicity rate calculated from catalog %s (%d events between t=[%.5lf, %.5lf]).** \n", catfile, sel_no, cat.tstart, cat.tend);
	}

	//find rate by counting events of required magnitude; is too few (<10), use GR instead.
	if (rate) *rate= cat.Z/T;		//fixme this does not take clustering into account. Should instead use value from below?
	if (minmag) *minmag= cat.Mc;


	if (crst.uniform){

		zlist=dvector(1,crst.nD);
		for (int i=1; i<=crst.nD; i++) zlist[i]=crst.depth[1+(i-1)*NP];

		rate_h=Helmstetter_cat(cat, crst, weights, ord);
		rate_v=fit_depth(zlist, zlist[2]-zlist[1], crst.nD, cat.depths0, cat.verr, weights, cat.Z);

		normv(rate_h, NP);
		normv(rate_v, crst.nD);

		*rate_grid=dvector(1,crst.N_allP);
		for (int p=1; p<=crst.N_allP; p++) {
			h_ind=(p-1)%NP+1;
			v_ind=(p-1)/NP+1;
			(*rate_grid)[p]=rate_h[h_ind]*rate_v[v_ind]*crst.N_allP;
		}
	}

	else{
		// todo [coverage] this block is never tested
		*rate_grid=Helmstetter_nonuni(crst.x, crst.y, crst.N_allP, cat.x0, cat.y0, cat.err, weights, cat.Z, ord);
	}

	free_cat(cat);
	if (zlist) free_dvector(zlist,1,crst.nD);
	if (rate_h)	free_dvector(rate_h, 1, NP);
	if (rate_v) free_dvector(rate_v, 1, crst.nD);

	return 0;

}

