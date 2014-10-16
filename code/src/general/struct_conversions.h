/*
 * readcatalog.h
 *
 *  Created on: Dec 22, 2011
 *      Author: camcat
 */


#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "../defines.h"
#include "../geom/coord_trafos.h"
#include "../util/nrutil.h"
#include "eqkfm_copy.h"
#include "mem_mgmt.h"


void eqk_filter(struct eqkfm **eqkfm1, int *Ntot, double Mag, double Depth);
int * cat_filter(struct catalog *cat, double Mag, double);
int *combine_eqkfm(struct eqkfm *eqkfm1, struct eqkfm *eqkfm2, int N1, int N2, double dt, double dM, double dR, int overwrite);
int *combine_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM);
double **union_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM, int***, int *);
double *timesfromeqkfm(struct eqkfm *eqkfm1, int N, int *);
double *magssfromeqkfm(struct eqkfm *eqkfm1, int N, int *);
double *timesfrompscmp(struct pscmp *DCFS, int N);
double *magsfrompscmp(struct pscmp *DCFS, int N);
void eqkfm2dist(struct eqkfm *eqkfm1, double *lats, double *lons, double *depths, int N, int Ntot, int);

