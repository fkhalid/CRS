/*
 * convert_geometry.h
 *
 *  Created on: Sep 13, 2013
 *      Author: camcat
 */

#ifndef CONVERT_GEOMETRY_H_
#define CONVERT_GEOMETRY_H_

#include <stdio.h>

#include "../defines.h"
#include "../util/nrutil.h"

int convert_geometry(struct crust crst, double *old_v, double **new_v, int sum, int increase_resolution);

#endif /* CONVERT_GEOMETRY_H_ */
