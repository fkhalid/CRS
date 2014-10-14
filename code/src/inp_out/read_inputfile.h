/*
 * read_inputfile.h
 *
 *  Created on: Nov 21, 2013
 *      Author: camcat
 */

#ifndef READ_INPUTFILE_H_
#define READ_INPUTFILE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../util/error.h"
#include "../util/nrutil.h"

int read_inputfile(char *input_fname, char *outname, char *reftime_str, char *crust_file, char *fore_template,
		char *catname, char ***focmeccat, char *background_rate_file, char *background_rate_cat, char *fixedmecfile, char *slipmodelfile, char *afterslipmodelfile,
		char *model_parameters_file, char *Logfile, struct tm *reftime,
		double *Tstart, double *Tend, long *seed, char *cmb_format, int *num_fm);

int read_slipformecfiles(char *inputfile, char ***listfiles, int *nfiles);
int read_listslipmodel(char *input_fname, struct tm reftime, struct slipmodels_list *allslipmodels, double res, int is_afterslip);

#endif /* READ_INPUTFILE_H_ */
