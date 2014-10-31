//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>

#include <stdio.h>

#include "../defines.h"
#define csep_forecast(...) csep_forecast_general(__VA_ARGS__,1)
#define csep_cmbmap(...) csep_forecast_general(__VA_ARGS__,0)

#define DAY2SEC 86400;
void csep_forecast_general(char *filename, struct crust crst, double *rates, int original_resolution,  int use_mags);
void write_csep_forecast(char *filename, double *lats, double *lons, double *deps, double dlat, double dlon, double ddep,
		double *mags, double dmag, double *rates, double *mag_fact, int NG, int Nmag);
