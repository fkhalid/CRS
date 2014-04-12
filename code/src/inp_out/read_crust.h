/*
 * read_crust.h
 *
 *  Created on: Jan 16, 2013
 *      Author: camcat
 */

#ifndef READ_CRUST_H_
#define READ_CRUST_H_


#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../okada/prestress.h"
#include "../util/error.h"
#include "../util/nrutil.h"
#include "read_csep_template.h"
#include "read_matrix.h"

int read_crust(char *fname, char *fnametemplate, char *focmecgridfile, struct crust *crst, double, double);
int read_farfalle_crust(char * file, struct crust *crst);
int read_pscmp_crust(char *fname, struct crust *crst);
int read_focmecgridfile(char *fname, struct crust *crst);

#endif /* READ_CRUST_H_ */

