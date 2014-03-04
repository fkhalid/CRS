/*
 * splines_eqkfm.h
 *
 *  Created on: May 16, 2013
 *      Author: camcat
 */

#ifndef SPLINES_EQKFM_H_
#define SPLINES_EQKFM_H_


#endif /* SPLINES_EQKFM_H_ */

#include "fit_splines.h"
#include "nrutil.h"
#include "../defines.h"
#include "../seis/soumod1.h"

void splines_eqkfm(struct eqkfm *eqkfm_aft, int Nas, int NF, double *times1, double *times2, int L, struct eqkfm *all_afterslips, long *seed);
