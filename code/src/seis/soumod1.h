/*
 * soumod-1.h
 *
 *  Created on: Mar 13, 2013
 *      Author: camcat
 */

#ifndef SOUMOD1_H_
#define SOUMOD1_H_

#endif /* SOUMOD1_H_ */

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../inp_out/print_output.h"
#include "../util/nrutil.h"
#include "../util/ran1.h"

int scale_to_mag(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double * slips, double *rakes);
int suomod1_taper(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, int top, int bottom, int right, int left);
int suomod1_resample(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double disc);

