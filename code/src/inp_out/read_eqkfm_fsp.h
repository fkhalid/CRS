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

int read_fsp_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF_out);
void track_position(long *pos, int NP, FILE* fin);
int next_separator(FILE * fin, char *string);
int find_key(FILE *fin, char *string, double *value);
int scan_nth(char *string, int n, double *result);
int read_slipvalues(FILE *fin, struct eqkfm *eqfm);


#endif /* READ_EQKFM_H_ */
