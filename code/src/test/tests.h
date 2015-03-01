/*
 * tests.h
 *
 *  Created on: Aug 23, 2013
 *      Author: camcat
 */

#ifndef TESTS_H_
#define TESTS_H_

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "../defines.h"
#include "../general/calculateDCFSperturbed.h"
#include "../general/forecast_stepG.h"
#include "../general/mem_mgmt.h"
#include "../general/setup.h"
#include "../geom/convert_geometry.h"
#include "../geom/coord_trafos.h"
#include "../geom/find_gridpoints.h"
#include "../inp_out/print_output.h"
#include "../inp_out/read_crust.h"
#include "../inp_out/read_csep_template.h"
#include "../inp_out/read_eqkfm.h"
#include "../inp_out/read_inputfile.h"
#include "../inp_out/read_matrix.h"
#include "../inp_out/read_zmap.h"
#include "../inp_out/write_csep_forecast.h"
#include "../okada/okadaDCFS.h"
#include "../okada/prestress.h"
#include "../seis/cmbopt.h"
#include "../seis/decluster.h"
#include "../seis/GR.h"
#include "../seis/Helmstetter.h"
#include "../seis/soumod1.h"
#include "../util/error.h"
#include "../util/hash.h"
#include "../util/moreutil.h"
#include "../util/nr.h"
#include "../util/nrutil.h"

int testspeed_coeff();
int test_hash();

#endif /* TESTS_H_ */
