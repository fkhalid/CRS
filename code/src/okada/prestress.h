/*
 * prestress.h
 *
 *  Created on: Apr 11, 2013
 *      Author: camcat
 */

#ifndef PRESTRESS_H_
#define PRESTRESS_H_

//#include <math.h>

//#include "../defines.h"
//#include "moreutil.h"
//#include "nrutil.h"

void prestress(double s1, double s2, double s3, double strike, double dip, double rake, double p,double f, double ***s);
double **prestress_eigen(double *s, double *str, double *dip);

#endif /* PRESTRESS_H_ */
