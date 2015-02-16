/*
 * forecast_stepG.h
 *
 *  Created on: Jul 20, 2012
 *      Author: camcat
 */
#ifndef FORECAST_STEPG_H_
#define FORECAST_STEPG_H_

#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "../defines.h"
#include "../util/nrutil.h"
#include "struct_conversions.h"

int rate_state_evolution(struct catalog cat, double *times, double **cmpdata, struct pscmp *DCFS, double tt0, double tt1, double Asig, double ta, int points[], double *NeX, double *NeT, double *Rate_end, int N, int NTS, int Neq, double *gamma_init, double *back_rate, double *R, int last);

#endif /* FORECAST_STEPG_H_ */
