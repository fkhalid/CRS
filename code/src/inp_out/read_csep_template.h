/*
 * read_txttemplate.h
 *
 *  Created on: Oct 4, 2013
 *      Author: camcat
 */

#ifndef READ_TXTTEMPLATE_H_
#define READ_TXTTEMPLATE_H_
#include <math.h>
#include <stddef.h>

#include "../defines.h"
#include "../geom/convert_geometry.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "read_matrix.h"

int read_fmindex(struct crust crst, char *fname, int **fm_index, int *no_zones);
int read_rate(struct crust crst, char *fname, double **bg_rate, double *minmag);
int read_csep_template(char *fname, int *no_magbins, int *nlat, int *nlon, int *ndep, int *ng, double *dlat, double *dlon, double *ddep, double *dmag,
		double **lats, double **lons, double **deps, double **rate, double *minlat, double *maxlat, double *minlon, double *maxlon, double *mindep, double *maxdep,
		double *minmag, double *maxmag, int *uni);
int read_csep_template_uni(char *fname, int *no_magbins, int *nlat, int *nlon, int *ndep, int *ng, double *dlat, double *dlon, double *ddep, double *dmag,
		double **lats, double **lons, double **deps, double **rate, double *minlat, double *maxlat, double *minlon, double *maxlon, double *mindep, double *maxdep,
		double *minmag, double *maxmag);

#endif /* READ_TXTTEMPLATE_H_ */
