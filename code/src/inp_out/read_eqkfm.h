/*
 * read_eqkfm.h
 *
 *  Created on: Jan 16, 2013
 *      Author: camcat
 */

#ifndef READ_EQKFM_H_
#define READ_EQKFM_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "../general/mem_mgmt.h"
#include "../general/setup.h"
#include "../general/struct_conversions.h"
#include "../seis/soumod1.h"
#include "../seis/WellsCoppersmith.h"
#include "../util/files.h"
#include "../util/nrutil.h"

int eqkfm_addslipmodels(struct eqkfm *eqfm1, struct slipmodels_list all_slipmodels,
						struct eqkfm **eqfm_comb, int N1,
						int *Ncomb, int **nfout, double dt, double dmag, double res,
						struct crust crst, struct flags flags);

int focmec2slipmodel(struct crust crst, struct eqkfm *eqfm1, double res, int refine, int taper);
int read_eqkfm(char *fname, char *cmbformat, struct eqkfm **eqfm1, int *NF_out, double *Mw, double mu);
int read_farfalle_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF_out);
int read_pscmp_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF2);

#endif /* READ_EQKFM_H_ */
