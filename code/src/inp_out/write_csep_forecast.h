//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>

#include <stdio.h>

#include "../defines.h"

#define DAY2SEC 86400;
void csep_forecast(char *filename, struct crust crst, double *rates);
void write_csep_forecast(char *filename, double *lats, double *lons, double *deps, double dlat, double dlon, double ddep,
		double *mags, double dmag, double *rates, double *mag_fact, int NG, int Nmag);
