/*
 * mscorr. *
 *  Created on: Apr 25, 2013
 *      Author: camcat
 *      translated from fortran code pscmp (Ronjiang Wang).
 *
 */

#include "mscorr.h"

double mscorr(double st1,double di1, double ra1, double st2, double di2, double ra2){

//corralation between two source mechanism (st,di and ra in degree)
      double ncorr,tcorr;
      double st[2],di[2],ra[2],ns[3][2],ts[3][2],rst[3],rdi[3];

      st[0]=st1*DEG2RAD;
      di[0]=di1*DEG2RAD;
      ra[0]=ra1*DEG2RAD;
      st[1]=st2*DEG2RAD;
      di[1]=di2*DEG2RAD;
      ra[1]=ra2*DEG2RAD;

      for (int j=0; j<2; j++){
        ns[1][j]=sin(di[j])*cos(st[j]+0.5*PI);
        ns[2][j]=sin(di[j])*sin(st[j]+0.5*PI);
        ns[3][j]=-cos(di[j]);
        rst[1]=cos(st[j]);
        rst[2]=sin(st[j]);
        rst[3]=0.0;
        rdi[1]=cos(di[j])*cos(st[j]+0.5*PI);
        rdi[2]=cos(di[j])*sin(st[j]+0.5*PI);
        rdi[3]=sin(di[j]);
        for (int i=0; i<3; i++) ts[i][j]=rst[i]*cos(ra[j])-rdi[i]*sin(ra[j]);
      }
      ncorr=0.0;
      tcorr=0.0;
      for (int i=0; i<3; i++){
        ncorr=ncorr+ns[i][1]*ns[i][2];
        tcorr=tcorr+ts[i][1]*ts[i][2];
      }
      return(ncorr*tcorr);
}
