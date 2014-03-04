/*
 * fit_splines.h
 *
 *  Created on: May 14, 2013
 *      Author: camcat
 */

#ifndef FIT_SPLINES_H_
#define FIT_SPLINES_H_


#endif /* FIT_SPLINES_H_ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil2.h"
#include "gasdev.h"
#include "interp_quad.h"
#include "spline.h"

void fit_splines(double *t, double *t2, int TS, int TS2, int N, double **slip_before, double *slip_before_err, double ***slip_after, int early_inter_mode, long *seed);
