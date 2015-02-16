/*
 * figures_CompPaper.c
 *
 *  Created on: Mar 27, 2014
 *      Author: camcat
 */

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../general/find_timesteps.h"
#include "../general/forecast_stepG.h"
#include "../general/mem_mgmt.h"
#include "../util/nrutil.h"

/* Compare numerical solution based on piecewise linear equation to analytical solution for logarithmic stress.
 * Analytical solution based on Dieterich (1994).
 */

void log_afterslip(char * Dtau_char){

	double t0=0.0, t1=300;
	double w=0.07;
	double Asig=10000;
	double ta=10000;
	double tau_coseismic=50000;
	double tau_postseismic=20000;
	double u=tau_postseismic/log(1+w*t1);
	double taudot=Asig/ta;
	double gamma0[2];
	gamma0[1]=1/taudot;
	double p=0.6, c=0.001;
	int L=100000;
	double times[L];
	char fname[200];
	FILE *fout;
	struct catalog cat;
	struct pscmp DCFS;
	double **DCFScont;
	double *rate;

	double Dtau;
	int afterslip=1;

	sscanf(Dtau_char,"%lf", &Dtau);

	DCFS.m=9.0;
	DCFS.nsel=1;
	DCFS.t=0.0;
	DCFS.which_pts=ivector(1,1);
	DCFS.cmb=dvector(1,1);
	DCFS.which_pts[1]=1;
	DCFS.cmb[1]=tau_coseismic;

	times[0]=-1e-14;

	sprintf(fname, "tests/log_afterslip2_%.0f.dat", Dtau);

	fout=fopen(fname,"w");
	findtimestepsomori(0.0, 1e-14, t1, t0, t1, tau_postseismic, Dtau, p, c, times+1, NULL, &L);
	printf("Dtau=%.1e\t%s -> %d steps\n", Dtau, fname, L);
	for (int i=1; i<=L; i++) fprintf(fout,"%.5e\t", times[i]);
	fprintf(fout,"\n");

	init_cat1(&cat, L);
	cat.mag[1]=9.0;
	for (int i=1; i<=L; i++) {
		cat.ngrid[i]=1;
		cat.ngridpoints[i][1]=1;
		cat.weights[i][1]=1.0;
		cat.t[i]=times[i];
	}

	rate=dvector(0,L);
	DCFScont=dmatrix(0,L, 1,1);
	for (int i=1; i<=L; i++) {
		DCFScont[i][1]= afterslip? u*(log(w*times[i]+1)-log(w*times[i-1]+1)) : 0;
		fprintf(fout,"%.5e\t", DCFScont[i][1]);
	}
	fprintf(fout,"\n");
//	for (int i=L; i<0; i--) DCFScont[i][1]-=DCFScont[i-1][1];

	rate_state_evolution(cat, times, DCFScont, &DCFS, t0, t1, Asig, ta, NULL, NULL, NULL, NULL, 1, L+1, 1, gamma0, NULL, rate, 0);
	for (int i=1; i<=L; i++) fprintf(fout,"%.5e\t", rate[i]);
	fprintf(fout,"\n");
	fclose(fout);


	return;


}
