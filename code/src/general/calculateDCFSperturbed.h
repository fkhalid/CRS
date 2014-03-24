/*
 * calculateDCFSperturbed.h
 *
 *  Created on: Mar 20, 2013
 *      Author: camcat
 */

#ifndef CALCULATEDCFSPERTURBED_H_
#define CALCULATEDCFSPERTURBED_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../okada/okadaDCFS.h"
#include "../seis/cmbopt.h"
#include "../seis/soumod1.h"
#include "../seis/WellsCoppersmith.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "../util/ran1.h"
#include "mem_mgmt.h"

void calculateDCFSperturbed(double **DCFSrand, struct pscmp *DCFS, struct eqkfm *eqkfmAf, struct eqkfm *eqkfm0, struct eqkfm *eqkfm1, struct flags flag,
		double *tevol, double *times, int Nmain, struct crust crst, struct Coeff_LinkList *AllCoeff, int NTScont, int NTSdisc, double **focmec,  int *fmzoneslim,
		int NFM, long *seed, int *fm_number, double tdata0, double tdata1, double Hurst, int refresh, int which_recfault);
void smoothen_DCFS(struct pscmp DCFS, int, int, int, long *seed, int, int **);
void smoothen_vector(int NgridT, int nLat, int nLon, int nD, double *values, long *seed, int**, int);

#endif /* CALCULATEDCFSPERTURBED_H_ */
