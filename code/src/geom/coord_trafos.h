/******************************************************************************
  "coord_trafos.c"

  description:  functions to transform geographical coordinates into the
                transversal mercator projection and vice versa


  project:  Simulation and parameter estimation using ETAS and shakemaps
   author:  Christoph Bach
     date:  2010-07-27
******************************************************************************/

#include <math.h>

/* transformation of geographical coordinates to local cartesian coordinates
 *
 * lat:  latitude
 * lon:  longitude
 * lat0: latitude of reference point
 * lon0: longitude of reference point
 * northern: north coordinate in [km] of new coordinate system
 * eastern: east coordinate in [km] of new coordinate system
 */
void latlon2localcartesian(double lat, double lon, double lat0, double lon0, double *northern, double *eastern);


/* transformation of local cartesian coordinates to geographical coordinates
 *
 * lat:  latitude
 * lon:  longitude
 * lat0: latitude of reference point
 * lon0: longitude of reference point
 * northern: north coordinate in [km] of new coordinate system
 * eastern:  east coordinate in [km] of new coordinate system
 */
void localcartesian2latlon(double northern, double eastern, double lat0, double lon0, double *lat, double *lon);
