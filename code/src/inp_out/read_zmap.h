/*
 * READ_ZMAP.H
 *
 *  Created on: Sep 24, 2013
 *      Author: camcat
 */

#ifndef READ_ZMAP_H_
#define READ_ZMAP_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../geom/find_gridpoints.h"
#include "../seis/GR.h"
#include "../util/nrutil.h"
#include "read_matrix.h"

int readZMAP (struct catalog *cat, struct eqkfm **eqfm, int *, char *file, struct crust crst, struct tm reftime, double t0s, double t1s,
		double t0c, double t1c, double Mmain, double tw, double border, double, double dDCFS, int findgridpoints);

int read_firstlineZMAP(char *file, struct tm reftime, double *time);

#endif /* READ_ZMAP_H_ */
