
/*   Copyright (C) 2015 by Camilla Cattania and Fahad Khalid.
 *
 *   This file is part of CRS.
 *
 *   CRS is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CRS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CRS.  If not, see <http://www.gnu.org/licenses/>.
 */


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

#include "../util/util1.h"

/* Compare numerical solution based on piecewise linear equation to analytical solution for logarithmic stress.
 * Analytical solution based on Dieterich (1994).
 */

void log_afterslip(){

	double t0=0.0, t1=300;
	double w=0.07;
	double Asig=100000;
	double ta=100000;
	double tau_coseismic=500000;
	double tau_postseismic=200000;
	double u=tau_postseismic/log(1+w*t1);
	double taudot=Asig/ta;
	double gamma0[2];
	double K;	//used for formula in paper.
	gamma0[1]=1/taudot;
	double p=0.6, c=0.001;
	int L=1000000;
	double times[L];
	char fname[200];
	FILE *fout, *fout0;
	struct catalog cat;
	struct pscmp DCFS;
	double **DCFScont;
	double *rate;

	int N1=3, N2=4;
	double Dtaus1[3]={1, 2, 5};
	double Dtaus2[4]={1, 2, 3, 4};
	double Dtau;
	int afterslip=1;

	DCFS.m=9.0;
	DCFS.nsel=1;
	DCFS.t=0.0;
	DCFS.which_pts=iarray(1,1);
	DCFS.cmb=darray(1,1);
	DCFS.which_pts[1]=1;
	DCFS.cmb[1]=tau_coseismic;

	times[0]=-1e-14;

	sprintf(fname, "tests/log_afterslip3_Ks.dat");
	fout0=fopen(fname,"w");

	for (int ii=0; ii<N1; ii++){
		for (int j=0; j<N2; j++){

			if (ii==0 && j==0) {
				Dtau=1454.5;//repeat for value used in code (K=0.177)
			}
			else{
				Dtau=Dtaus1[ii]*pow(10.0, Dtaus2[j]);
			}

			L=1000000;

			sprintf(fname, "tests/log_afterslip3_%.0f.dat", Dtau);

			fout=fopen(fname,"w");
			//findtimestepsomori now takes K as input.
			//findtimestepsomori(0.0, 1e-14, t1, t0, t1, tau_postseismic, Dtau, p, c, times+1, &K, &L);

			fprintf(fout0, "%.3lf\t %.3lf\n", Dtau, K);	//print out to summary file.

			printf("Dtau=%.1e\t%s -> %d steps\n", Dtau, fname, L);
			for (int i=1; i<=L; i++) fprintf(fout,"%.5e\t", times[i]);
			fprintf(fout,"\n");

			init_cat1(&cat, L);
			cat.mag[1]=9.0;
			for (int i=1; i<=L; i++) {
				cat.ngridpoints[i]=malloc(2*sizeof(double));
				cat.weights[i]=malloc(2*sizeof(double));
				cat.ngrid[i]=1;
				cat.ngridpoints[i][1]=1;
				cat.weights[i][1]=1.0;
				cat.t[i]=times[i];
			}

			rate=darray(0,L);
			DCFScont=d2array(0,L, 1,1);
			for (int i=1; i<=L; i++) {
				DCFScont[i][1]= afterslip? u*(log(w*times[i]+1)-log(w*times[i-1]+1)) : 0;
				fprintf(fout,"%.5e\t", DCFScont[i][1]);
			}
			fprintf(fout,"\n");
		//	for (int i=L; i<0; i--) DCFScont[i][1]-=DCFScont[i-1][1];

			//rate_state_evolution(cat, times, DCFScont, &DCFS, t0, t1, t1-t0, Asig, ta, NULL, NULL, NULL, NULL, 1, L+1, 1, gamma0, NULL, rate, 0);
			for (int i=1; i<=L; i++) fprintf(fout,"%.5e\t", rate[i]);
			fprintf(fout,"\n");
			fclose(fout);
		}
	}

	fclose(fout0);

	return;


}
