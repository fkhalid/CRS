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
#include "../geom/dft2d.h"
#include "../inp_out/print_output.h"
#include "../util/nrutil.h"
#include "../util/ran1.h"

int scale_to_mag(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double * slips, double *rakes);
int adjust_faults(struct eqkfm *eqkfm0,  int NF, int vert);
int which_taper(struct eqkfm *eqkfm0,  int NF, int tap_bot, int tap_top, double);
int suomod1_taper(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
int suomod1_resample(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double disc, double velmean);
int suomod1_hf(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double velmean, long *seed, int);
int suomod1_add(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double velmean, long *seed);
int suomod1_addnoise(struct eqkfm eqkfm1, struct eqkfm eqkfm2, struct eqkfm *eqkfm);
int suomod1_cleanup(struct eqkfm *eqkfm2);
//int suomod1_addhf(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double velmean, long *seed);

