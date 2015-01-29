/*
 * background_rate.h
 *
 *  Created on: Nov 12, 2013
 *      Author: camcat
 */

#ifndef BACKGROUND_RATE_H_
#define BACKGROUND_RATE_H_

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../inp_out/read_matrix.h"
#include "../inp_out/read_zmap.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "decluster.h"
#include "GR.h"
#include "Helmstetter.h"

int background_rate(char *catfile, struct crust *crst_in, struct tm reftime,
		double Mmain, double *minmag, double *rate, double **rate_grid, double dR, double dZ, double min_smoothing, int ord);

#endif /* BACKGROUND_RATE_H_ */
