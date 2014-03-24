//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>

//#include "../defines.h"
//#include "nrutil.h"

#define DAY2SEC 86400;

#include <stdio.h>
#include <time.h>

#include "../defines.h"
#include "../util/nrutil.h"

//void printCSEPforecast(char *filename, double *lats, double *lons, double *deps, double *mags, double *rates, double *mag_fact, int NG, int Nmag);
int sum_DCFS(struct pscmp *DCFS, double **cmb, int N, int Ntot);
int print_rate(char *fname, struct crust crst, double Mc, double *rate);
int print_grid(char *fname, struct pscmp DCFS, struct crust, double *rate);
int print_slipmodel(char* filename, struct eqkfm *eqfm1, int NF);
int print_cat(char *fname, struct catalog cat);
//void printXMLforecast(char *filename, double *lats, double *lons, double *deps,  int Nlat, int Nlon, int Ndep, double Dmax, double Dmin, double *mags,
//		double *rates, double *mag_fact, int Nmag, double t0, double t1, struct tm nowtime, struct tm eqktime);
