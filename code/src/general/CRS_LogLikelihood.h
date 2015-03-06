/*
 * CRS_LogLihelihood.h
 *
 *  Created on: Jul 23, 2013
 *      Author: camcat
 */

#ifndef CRS_LOGLIHELIHOOD_H_
#define CRS_LOGLIHELIHOOD_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../geom/convert_geometry.h"
#include "../inp_out/print_output.h"
#include "../inp_out/write_csep_forecast.h"
#include "../util/nrutil.h"
#include "calculateDCFSperturbed.h"
#include "forecast_stepG.h"

int CRSforecast (double *LL, int Nsur, struct pscmp *DCFS, struct eqkfm *eqkfm_aft, struct eqkfm *eqkfm0, struct flags flags,
		struct crust crst, struct Coeff_LinkList *AllCoeff, int NTScont, int Nm, int Na, int NgridT, double **focmec, int *fmzonelim, int NFM,
		long *seed, struct catalog cat, double *times, double tstart, double t0, double t1, double dtstep, double Asig, double ta, double r0,
		double **all_gammas0, int multiple_input_gammas, int fromstart,
		char * print_cmb, char *print_forex, char *print_foret, char * printall_cmb, char *printall_forex, char *printall_foret, char *print_LL);

int CRSLogLikelihood (double *LL, double *Ldum0_out, double *Nev, double *I, double *r_out, int Nsur, struct pscmp *DCFS,
		struct eqkfm *eqkfm_aft, struct eqkfm *eqkfm0, struct flags flags,
		struct crust crst, struct Coeff_LinkList *AllCoeff, int NTScont, int Nm, int Na, int NgridT, double **focmec, int *fmzonelim, int NFM,
		long *seed, struct catalog cat, double *times, double tstart, double tt0, double tt1, double tw, double Mag_main, double Asig, double ta, double r0, int fixr,
		double *gammas0, double **all_new_gammas, int fromstart,
		char * printall_cmb, char *printall_forex, int refresh);

#ifdef _CRS_MPI
// [Fahad] -- Added for use in partition size calculation when using
//			  MPI enabled code.
inline int roundUpFrac(double x) {
	double integralPart, fractionalPart, result;

	fractionalPart = modf(x, &integralPart);

	if(integralPart != 0.0 && fractionalPart != 0.0) {
		result = ceil(integralPart + 0.5);
	}
	else {
		result = integralPart;
	}

	return result;
}
#endif


#endif /* CRS_LOGLIHELIHOOD_H_ */
