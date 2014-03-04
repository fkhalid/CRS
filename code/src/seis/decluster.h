/*
 * decluster.h
 *
 *  Created on: Nov 11, 2013
 *      Author: camcat
 */

#ifndef DECLUSTER_H_
#define DECLUSTER_H_

#include <math.h>

#include "../defines.h"
#include "../util/nrutil.h"

int *decluster_catalog_rescalegrid(struct catalog cat, struct crust crst, double Mmain, double **time_missing, int d3);
int * decluster_catalog(struct catalog cat, double Mmain, double **time_missing, int d3);
void KG74(double M, double *D, double *T);

#endif /* DECLUSTER_H_ */
