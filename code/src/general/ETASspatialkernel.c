/*
 * ETASspatialkernel.c
 *
 *  Created on: Jun 19, 2014
 *      Author: camcat
 */

#include <math.h>

#include "ETASspatialkernel.h"
#include "../defines.h"
#include "struct_conversions.h"

void ETASspatialkernel(struct pscmp *DCFS, struct eqkfm *eqkfm1, struct crust crst, int NTSdisc, double tdata0, double tdata1,
		double etas_d, double etas_q, double etas_alpha, double Mc){

	eqkfm2dist(eqkfm1, crst.lat, crst.lon, crst.depth, crst.N_allP, NTSdisc, 1);	//may be duplicate but just to be sure (and don't want to change previous non-ETAS code).

	for (int eq1=0; eq1<NTSdisc; eq1++){
		if (DCFS[eq1].t <tdata0 || DCFS[eq1].t>tdata1) continue;
		if (DCFS[eq1].nsel!=0){
			isoETAS(DCFS[eq1], etas_d, etas_q, etas_alpha, Mc);
		}
	}
}

//--------------------------------------------------------------//

int isoETAS(struct pscmp DCFS, double d, double q, double alpha, double Mc){
	double r;
	int Nsel;
	double tot=0.0;

	Nsel=DCFS.nsel;
	double totno=pow(10.0,(DCFS.m-Mc)*alpha);
	#pragma omp parallel for private(r) reduction(+: tot)
	for (int i=1; i<=Nsel; i++){
		r=DCFS.fdist[i];	//r is in km.
		//DCFS.cmb[i]=totno*fspatial(r, q, d);
		DCFS.cmb[i]=fspatial(r, q, d);
		tot+=DCFS.cmb[i];
	}
	for (int i=1; i<=Nsel; i++) DCFS.cmb[i]*=totno/tot;

	return(0);
}

double fspatial(double rr, double q, double d) {

	double res;
	double a6m1 = q - 1.0;
	double x2a6m1 = 2.0 * a6m1;
	double sqa7 = sq(d);
	double sqrr = sq(rr);
	double sqrrsqa7 = sqrr + sqa7;
	double pow1 = pow(sqrrsqa7, q);
	double pow2 = pow(d, x2a6m1);
	double pow3 = pow2 / pow1;

	res = a6m1 / PI * pow3; // (q-1)/pi * d^(2(q-1))/(r^2+d^2)^q

	return res;
}
