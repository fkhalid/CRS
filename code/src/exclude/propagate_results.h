/*
 * propagate_results.h
 *
 *  Created on: Aug 30, 2013
 *      Author: camcat
 */

#ifndef PROPAGATE_RESULTS_H_
#define PROPAGATE_RESULTS_H_

#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../util/nrutil.h"
#include "read_matrix.h"

int check_if_snapshot_exists(char * folder, double *t, struct tm reftime, long *);
int load_oldLL(char * folder, double ***LLs);
int check_if_no_newdata();
void write_gammas(char *folder, int p, double **gammas, int Nsur, int NG);
void load_gammas(char * folder, int p, double **gammas_old, int NG);

#endif /* PROPAGATE_RESULTS_H_ */

