/*
 * lin_interp_eqkfm.h
 *
 *  Created on: May 29, 2015
 *      Author: camcat
 */

#ifndef LIN_INTERP_EQKFM_H_
#define LIN_INTERP_EQKFM_H_

void fit_lin(double *times1, double *times2, int Nas, int L, int *ind, int NP, double **slipbefore_st, double **slip_after);
void lin_interp_eqkfm(struct eqkfm **eqkfm_aft, int NF, int L, double *times2, int *ind);

#endif /* LIN_INTERP_EQKFM_H_ */
