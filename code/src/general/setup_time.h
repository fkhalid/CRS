/*
 * setup_time.h
 *
 *  Created on: May 29, 2015
 *      Author: camcat
 */

#ifndef SETUP_TIME_H_
#define SETUP_TIME_H_

#include "../defines.h"

int timesteps_log(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2, double smallstepstime, double TAU, double dtau, double timeTAU);

int timesteps_lin(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2, int ***allind);

int setup_aseismic_multi_linear(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2);

int setup_aseismic_single_linear(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2);

int setup_aseismic_splines(double t0, double t1, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2, long *seed);

int setup_aseismic_single_log(double t0, double t1, double ts, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2, long *seed);

#endif /* SETUP_TIME_H_ */
