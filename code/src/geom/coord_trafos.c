/******************************************************************************
  "coord_trafos.c"

  description:  functions to transform geographical coordinates into the
                transversal mercator projection and vice versa


  project:  Simulation and parameter estimation using ETAS and shakemaps
   author:  Christoph Bach
     date:  2010-07-27
******************************************************************************/

#include "coord_trafos.h"

void latlon2localcartesian(double lat, double lon, double lat0, double lon0,
                           double *northern, double *eastern) {

  double e, e2, N, T, C, A, M, M0;
  double phi, lambda, phi0, lambda0, dlambda;

  double k0 = 0.9996;      /* scale on central meridian for UTM */
  double a = 6378137.0;    /* equatorial radius */
  double b = 6356752.3142; /* polar radius */
  double pi = 3.141592653589793;

  /* phi and lambda have to be in radian */
  phi     = lat * pi/180.0;
  lambda  = lon * pi/180.0;
  phi0    = lat0 * pi/180.0;
  lambda0 = lon0 * pi/180.0;

  /* eccentricity of the earth's elliptical cross-section */
  e = sqrt(((a * a) - (b * b)) / (a * a));
  /* e'² */
  e2 = pow(e, 2.0) / (1.0 - pow(e, 2.0));

  N = a / sqrt(1.0 - pow(e, 2.0) * pow(sin(phi), 2.0));
  T = tan(phi);
  C = e2 * pow(cos(phi), 2.0);
  dlambda=lambda - lambda0;
  if (fabs(dlambda)>pi) dlambda= (dlambda>0)? dlambda-2*pi : dlambda+2*pi;
  A = dlambda * cos(phi);
  /* M is the true distance along the central meridian from the Equator to phi */
  M  = a*((1.0 - pow(e,2.0)/4.0 - 3.0*pow(e,4.0)/64.0 - 5.0*pow(e,6.0)/256.0)*phi
      - (3.0*pow(e,2.0)/8.0 + 3.0*pow(e,4.0)/32.0 + 45.0*pow(e,6.0)/1024.0)*sin(2.0*phi)
      + (15.0*pow(e,4.0)/256.0 + 45.0*pow(e,6.0)/1024.0)*sin(4.0*phi)
      - (35.0*pow(e,6.0)/3072.0)*sin(6.0*phi));
  /* M0 = M calculated for phi0 */
  M0 = a*((1.0 - pow(e,2.0)/4.0 - 3.0*pow(e,4.0)/64.0 - 5.0*pow(e,6.0)/256.0)*phi0
      - (3.0*pow(e,2.0)/8.0 + 3.0*pow(e,4.0)/32.0 + 45.0*pow(e,6.0)/1024.0)*sin(2.0*phi0)
      + (15.0*pow(e,4.0)/256.0 + 45.0*pow(e,6.0)/1024.0)*sin(4.0*phi0)
      - (35.0*pow(e,6.0)/3072.0)*sin(6.0*phi0));

  /* new coordinates in [km] relative to lat0/lon0 */
  *northern = (k0*(M - M0 + N*tan(phi)*(pow(A,2)/2 + (5 - T + 9*C + 4*pow(C,2))*pow(A,4)/24
              + (61 - 58*T + pow(T,2) + 600*C - 330*e2)*pow(A,6)/720)))/1000;
  *eastern  = (k0*N*(A + (1 - T + C)*pow(A,3)/6
              + (5 - 18*T + pow(T,2) + 72*C - 58*e2)*pow(A,5)/120))/1000;

}

//void localcartesian2latlon(double northern, double eastern, double lat0, double lon0, double *lat, double *lon){
//  double M0, e, e2, M, mu, e1, phi1, C1, T1, N1, R1, D, phi0, lambda0, phi, lambda;
//
//  double k0 = 0.9996;        /* scale on central meridian for UTM */
//  double a = 6378137.0;      /* equatorial radius */
//  double b = 6356752.3142;   /* polar radius */
//  double pi = 3.141592653589793;
//
//  phi0    = lat0 * pi/180.0;
//  lambda0 = lon0 * pi/180.0;
//
//  /* convert [km] to [m] */
//  eastern *= 1000.0;
//  northern *= 1000.0;
//
//  /* eccentricity of the earth's elliptical cross-section */
//  e = sqrt(((a * a) - (b * b)) / (a * a));
//  /* e'² */
//  e2 = pow(e, 2.0) / (1.0 - pow(e, 2.0));
//
//
//  /* M0 = M calculated for phi0 */
//  M0 = a*((1 - pow(e,2.0)/4.0 - 3.0*pow(e,4.0)/64.0 - 5.0*pow(e,6.0)/256.0)*phi0
//      - (3.0*pow(e,2.0)/8.0 + 3.0*pow(e,4.0)/32.0 + 45*pow(e,6.0)/1024.0)*sin(2.0*phi0)
//      + (15.0*pow(e,4.0)/256.0 + 45.0*pow(e,6.0)/1024.0)*sin(4.0*phi0)
//      - (35.0*pow(e,6.0)/3072.0)*sin(6.0*phi0));
//
//  M = M0 + northern/k0;
//
//  mu = M/(a*(1 - pow(e,2)/4 - 3*pow(e,4)/64 - 5*pow(e,6)/256));
//  e1 = (1 - sqrt(1 - pow(e,2)))/(1 + sqrt(1 - pow(e,2)));
//  phi1 = mu + (3*e1/2 - 27*pow(e1,3)/32)*sin(2*mu) + (21*pow(e1,2)/16 - 55*pow(e1,4)/32)*sin(4*mu) + (151*pow(e1,3)/96)*sin(6*mu) + (1097*pow(e1,4)/512)*sin(8*mu);
//
//  C1 = e2*pow(cos(phi1),2.0);
//  T1 = pow(tan(phi1),2.0);
//  N1 = a/sqrt(1.0 - pow(e,2)*pow(sin(phi1),2.0));
//  R1 = a*(1.0 - pow(e,2.0))/pow((1.0 - pow(e,2)*pow(sin(phi1),2.0)),(3.0/2.0));
//  D  = eastern/(N1*k0);
//
//
//  phi = phi1 - (N1*tan(phi1)/R1)*(pow(D,2.0)/2.0 - (5.0 + 3.0*T1 + 10.0*C1 - 4.0*pow(C1,2.0) - 9.0* e2)*pow(D,4.0)/24.0 + (61.0 + 90.0*T1 + 298.0*C1 + 45.0*pow(T1,2.0) - 252.0*e2 - 3.0*pow(C1,2.0))*pow(D,6.0)/720.0);
//  lambda = lambda0 + (D - (1 + 2.0*T1 + C1)*pow(D,3.0)/6.0 + (5.0 - 2.0*C1 + 28.0*T1 - 3.0*pow(C1,2.0) + 8.0*e2 + 24.0*pow(T1,2.0))*pow(D,5.0)/120.0)/cos(phi1);
//
//  *lat = phi * 180.0/pi;
//  *lon = lambda * 180.0/pi;
//
//}

