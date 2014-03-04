/*
 * dft2d.h
 *
 *  Created on: Mar 12, 2013
 *      Author: camcat
 */

#ifndef DFT2D_H_
#define DFT2D_H_

#include <math.h>

#include "../util/nrutil.h"

void dft2d(int n, int m, int inverse, double *gRe, double *gIm, double *GRe, double *GIm);

#endif /* DFT2D_H_ */

