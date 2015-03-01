/*
 * ETASspatialkernel.h
 *
 *  Created on: Jun 19, 2014
 *      Author: camcat
 */

#ifndef ETASSPATIALKERNEL_H_
#define ETASSPATIALKERNEL_H_

#endif /* ETASSPATIALKERNEL_H_ */
/*
 * ETASspatialkernel.c
 *
 *  Created on: Jun 19, 2014
 *      Author: camcat
 */

#include <math.h>
#include "../defines.h"

void ETASspatialkernel(struct pscmp *DCFS, struct eqkfm *eqkfm1, struct crust crst, int NTSdisc, double tdata0, double tdata1, double etas_d, double etas_q, double etas_alpha, double Mc);
int isoETAS(struct pscmp DCFS, double d, double q, double alpha, double Mc);
double fspatial(double rr, double q, double d);
