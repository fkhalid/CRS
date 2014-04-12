/*
 * Helmstetter.h
 *
 *  Created on: Nov 8, 2013
 *      Author: camcat
 */

#ifndef HELMSTETTER_H_
#define HELMSTETTER_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../geom/find_gridpoints.h"
#include "../util/nrutil.h"

double * Helmstetter(double *xgrid, double *ygrid, double dx, double dy, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord);
double *Helmstetter_nonuni(double *xgrid, double *ygrid, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord);
double *Helmstetter_cat(struct catalog cat, struct crust crst, double *weights, int ord);
double *fit_depth(double *zgrid, double dz, int Ngrid, double *zs, double *err, double *weights, int N);

#endif /* HELMSTETTER_H_ */
