/*
 * splines_eqkfm.h
 *
 *  Created on: May 16, 2013
 *      Author: camcat
 */

#ifndef SPLINES_EQKFM_H_
#define SPLINES_EQKFM_H_

void splines_eqkfm(struct eqkfm **eqkfm_aft, int Nas, int NF, double *times1, double *times2, int L, long *seed);
//void splines_eqkfm(struct eqkfm *eqkfm_aft, int Nas, int NF, double *times1, double *times2, int L, struct eqkfm *all_afterslips, long *seed);

#endif /* SPLINES_EQKFM_H_ */

