/*
 * interp_quad.c
 *
 *  Created on: May 14, 2013
 *      Author: camcat
 */

#include "interp_quad.h"

void interp_quad(double t0, double t1, double y, double y1, double *dts, double *values, int ndts){
	//t1= time to which y, y1 refer; t0=time at which y=0;
	//dts are time steps at which function is evaluated (ndts of them).
	// y, y1 contain the y value and first derivative at t=t1.

	double DT=t1-t0;
	double a,b;	//polinomial coefficients.

	a=(y1*DT-y)/pow(DT,2);
	b=2*y/DT - y1;
	for (int v=1; v<=ndts; v++) values[v]= a*((dts[v]-t0)*(dts[v]-t0))+b*(dts[v]-t0);

	return;
}
