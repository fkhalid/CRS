/*
 * dist2fault.h
 *
 *  Created on: Apr 18, 2013
 *      Author: camcat
 */

#ifndef DIST2FAULT_H_
#define DIST2FAULT_H_


#include <math.h>

#include "../defines.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"

void dist2faultDCFS(struct pscmp *DCFS, struct crust, struct eqkfm *eqkfm1);
double *dist2fault0(double *lats, double *lons, double *depths, int NP, double strike, double dip, double Lat0, double Lon0, double D0, double *pos_s, double *pos_d);
void dist2fault(double *lats, double *lons, double *depths, double *dist, int NP, int NF, struct eqkfm *eqkfm1);
void dist2line(double *P, double *A, double *B, int d3, int finitesegment, double *D, double *d);
void select_onofffault(double *FMdist, int NFM, int **labove, int **lbelow, int **lon, int *, int *, int*);

#endif /* DIST2FAULT_H_ */
