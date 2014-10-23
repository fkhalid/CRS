/*
 * setup_eqkfm.h
 *
 *  Created on: Jul 24, 2013
 *      Author: camcat
 */

#ifndef SETUP_EQKFM_H_
#define SETUP_EQKFM_H_


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../defines.h"
#include "../geom/coord_trafos.h"
#include "../inp_out/read_eqkfm.h"
#include "../inp_out/read_focmec.h"
#include "../inp_out/read_matrix.h"
#include "../inp_out/read_zmap.h"
#include "../okada/okadaDCFS.h"
#include "../seis/soumod1.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "../util/splines_eqkfm.h"
#include "eqkfm_copy.h"
#include "find_timesteps.h"
#include "mem_mgmt.h"
#include "struct_conversions.h"

void set_current_slip_model(struct eqkfm *eqkfm0, int slipmodel_index);
int setup_catalogetc(char *catname, char **focmeccat, int nofmcat, struct tm reftime, double dDCFS, double Mag_main, struct crust crst,
		struct catalog *cat, struct eqkfm **eqkfm1, double ***focmec, int **firstelements, struct flags flag, int *NFM, int *Ntot, int *Nmain,
		double dt, double dM, double xytoll, double ztoll, double dR, double tw, double tstart, double tend);
int setup_afterslip_eqkfm(struct slipmodels_list slipmodels, struct crust crst, struct eqkfm **eqkfm0res);
int setup_eqkfm_element(struct eqkfm *eqkfm0res, char **slipmodel, char *cmb_format, int no_slipmodels, double mu, double disc,
		double tmain, int nsel, int *sel_pts, double *mmain, int cuts_surf, int *NF0, double lat0, double lon0);
int load_newdata(double *t0, double t1, struct set_of_models *allmodels, int Nmain, int *NFaults, char **slipmodels, char **multimodels, int *no_slipmodels, int *Nmain_now);
int setup_afterslip_evol(double Teq, double t0, double t1, double *Cs, double *ts, int Nfun, struct eqkfm **eq_aft, double *t_afterslip, int Nas,int Nfaults,
		int afterslip, int *L, double **times2, double **tevol_afterslip, long *seed);
int mask_afterslip(double time, double *times, int L, double *evol_afterslip0, double **evol_afterslip, int *Lnow);
int setup_CoeffsDCFS(struct Coeff_LinkList **Coefficients, struct pscmp **DCFS_out, struct crust crst, struct eqkfm *eqkfm0,
		int Nm, int *Nfaults);
int update_CoeffsDCFS(struct Coeff_LinkList **Coefficients, struct crust crst, struct eqkfm *eqkfm0, int Nm, int *Nfaults);
#endif /* SETUP_EQKFM_H_ */
