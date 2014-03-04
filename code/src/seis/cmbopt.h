/*
 * cmbopt.h
 *
 *  Created on: Apr 25, 2013
 *      Author: camcat
 */

#ifndef CMBOPT_H_
#define CMBOPT_H_


#endif /* CMBOPT_H_ */

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../okada/okadaDCFS.h"
#include "../util/mscorr.h"
#include "../util/nrutil.h"
#include "../util/roots3.h"

void cmbopt(double sxx, double syy, double szz, double sxy, double syz, double szx, double p, double f, double st0, double di0, double ra0, double *cmb,double *st1, double *di1, double *ra1, double *st2, double *di2, double *ra2);
void DCFScmbopt(struct pscmp *DCFS, int, struct crust crst);
