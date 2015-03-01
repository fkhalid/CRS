/*
 * ETASforecast.h
 *
 *  Created on: Jun 19, 2014
 *      Author: camcat
 */

#ifndef ETASFORECAST_H_
#define ETASFORECAST_H_

/*
 * ETASforecast.c
 *
 *  Created on: May 17, 2012
 *      Author: camcat
 */

#include "ETASforecast.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>

#include "../defines.h"
#include "../util/nrutil.h"
#include "struct_conversions.h"

#define Max 100000
#define MaxTS 1000
#define Maxeq 10000

int ETASforecast (struct pscmp *DCFS, struct eqkfm *eqkfm, struct crust crst, int NTSdisc, int NgridT, int Nm, struct catalog cat,
		double *tts, int Ntts, double tw, char *print_forex0, char *print_foret,  char *print_LL, double etas_p,  double etas_c,  double etas_K,  double etas_d,
		double etas_q,  double etas_mu, double etas_alpha, double Mc);

int forecastETAS(struct catalog cat, struct pscmp *DCFS, double tt0, double tt1, int points[], double *out_NeX, double *NeT,
		double *Rate_end, int N, int Neqks, double *R, double etas_c, double etas_p, double etas_K);

#endif /* ETASFORECAST_H_ */
