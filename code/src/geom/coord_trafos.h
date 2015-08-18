/******************************************************************************
  "coord_trafos.c"

  description:  functions to transform geographical coordinates into the
                transversal mercator projection and vice versa


  project:  Simulation and parameter estimation using ETAS and shakemaps
   author:  Christoph Bach
     date:  2010-07-27
******************************************************************************/

#include <math.h>

void latlon2localcartesian(double lat, double lon, double lat0, double lon0, double *northern, double *eastern);
