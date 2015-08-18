/*
 * findtimesteps.h
 *
 *  Created on: Mar 27, 2012
 *      Author: camcat
 */

#ifndef FINDTIMESTEPS_H_
#define FINDTIMESTEPS_H_

#include <math.h>
#include <stdio.h>
#include "../defines.h"

#include "../util/nrutil.h"

int findtimestepsomori(double te, double t0,double t1, double K, double p, double c, double *times, double *K_over_tau0, int *L);

#endif /* FINDTIMESTEPS_H_ */
