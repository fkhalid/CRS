/*
 * GR.h
 *
 *  Created on: Sep 20, 2013
 *      Author: camcat
 */

#ifndef GR_H_
#define GR_H_

#include <stdlib.h>
#include <math.h>
#include "../util/nrutil.h"

double *assign_GRnorm(double *mags, int N, double b, int Minf);
int compare (const void * a, const void * b);
int bin_equnumber(double *v, int N, int Nbin, double **bin_c, double **norm_count);
double Mc_maxcurv(double *mags, int N);
double calculatebvalue(double *mags, int N, double Mc);

#endif /* GR_H_ */
