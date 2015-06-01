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

void findtimestepslog(double t0,double t1, double tau0, double Dtau, double w, double **times,int *L);
int findtimestepsomori(double te, double t0,double t1, double K, double p, double c, double *times, double *K_over_tau0, int *L);
void tevolomori(double te, double *times, double Kotau, double p, double c, double *tevol, int L);

#endif /* FINDTIMESTEPS_H_ */
