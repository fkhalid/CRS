/*
 * readfocmec.h
 *
 *  Created on: Feb 17, 2012
 *      Author: camcat
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../geom/find_gridpoints.h"
#include "../seis/WellsCoppersmith.h"
#include "../util/nrutil.h"
#include "../general/eqkfm_copy.h"
#include "read_matrix.h"

int readmultiplefocmec(char **focmecfiles, int nofiles, char *which_format,
					   struct crust crst, double, double, double dDCFS,
					   struct tm reftime, double t0, double t1, double tfocmec,
					   double mag, double ***focmec, int **firstelements,
					   int *NFM, int *NFM_timesel, struct eqkfm **eqkfm,
					   int sel, int fm2);

int readfocmec(char *focmecfile, char *which_format, struct crust crst, double, double, double dDCFS, struct tm reftime,
		double t0, double t1, double tfocmec, double mag, double ***focmec,	int *NFM, int *NFM_timesel, struct eqkfm **eqkfm,int sel, int fm2);

void select_fm_time(double **focmec, int *NFM, double Tstart);
