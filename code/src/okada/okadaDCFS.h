/*
 * okadaDCFS.h
 *
 *  Created on: Jan 4, 2013
 *      Author: camcat
 */

#ifndef OKADADCFS_H_
#define OKADADCFS_H_

#include <math.h>
#include <omp.h>
#include <stdio.h>

#include "../defines.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "pscokada.h"

//----------------top level functions----------------//
int resolve_DCFS(struct pscmp DCFS, struct crust crst, double *strikeRs, double *dipRs, double *rake, int optrake);
int okadaDCFS(struct pscmp DCFS, struct eqkfm *eqkfm1, int NF, struct crust crst, double *strikeR, double *dipR, int full_tensor);
int okadaCoeff(float ****Coeffs_st, float ****Coeffs_dip, struct eqkfm *eqkfm1, int NF, struct crust crst, double *lats, double *lons, double *depths);
int okadaCoeff_mpi(float ****Coeffs_st, float ****Coeffs_dip, struct eqkfm *eqkfm1, int NF, struct crust crst, double *lats, double *lons, double *depths);
int okadaCoeff2DCFS(float ***Coeffs_st, float ***Coeffs_d, struct pscmp DCFS, struct eqkfm *eqkfm1, struct crust crst, double *strikeR, double *dipR, int full_tensor);
int okadaCoeff_resolve(struct Coeff_LinkList Coeffs, float ***Coeffs_st2, float ***Coeffs_di2, struct crust crst, double *strikeRs, double *dipRs, double *rakeRs);
int resolvedCoeff2DCFS(float **Coeffs_st, float **Coeffs_d, struct pscmp DCFS, struct eqkfm *eqkfm1, struct crust crst);
int isoDCFS(struct pscmp DCFS, struct eqkfm eqkfm1);

//----------------auxiliary functions----------------//
int choose_focmec(struct eqkfm eqkfm1, double *strike, double *dip, double *rake);
void patch_pos(struct eqkfm eqfm, int p, double *east, double *north, double *depth);
double ** comp2tensor(float *v, double ***S0);
double *normal_vector(double, double);
double *slip_vector(double strikeR, double dipR, double rakeR);
double *opt_s(double *stress, double sigma, double *n, double *result);
double *sum_v(double *v1, double *v2, double *sum, int N);
double resolve_S(double **S, double strikeR, double dipR, double rakeR, double f, double *stress0, double sigma0, double *newrake, int opt_rake);
double resolve_n(double **S, double *n, double *rake, double fric, double *stress0, double sigma0, double *slip_v);

#endif /* OKADADCFS_H_ */
