//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>

#include <stdio.h>

#include "../defines.h"

// Define macros to print out forecast (no. of earthquakes, binned into magnitudes) and cmbmap (stress field, not binned into magnitudes).
#define csep_forecast(...) csep_forecast_general(__VA_ARGS__,1)
#define csep_cmbmap(...) csep_forecast_general(__VA_ARGS__,0)

void csep_forecast_general(char *filename, struct crust crst, double *rates, int original_resolution,  int use_mags);
void write_csep_forecast(char *filename, double *lats, double *lons, double *deps, double dlat, double dlon, double ddep,
		double *mags, double dmag, double *rates, double *mag_fact, int NG, int Nmag);
