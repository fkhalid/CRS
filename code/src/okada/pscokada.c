#include "pscokada.h"

#include <math.h>
//#include <stdio.h>

#include "../defines.h"
#include "../util/nrutil.h"
#include "dc3d.h"



/************************************************************/
/*			cmbfix				    */
/*		calculate Coulomb Stress		    */
/*							    */
/* input:						    */
/*  stress tensor, pore pressure, friction coefficient	    */
/*  receiver orientation parameter (strike, dip and rake)   */
/*							    */
/* return:						    */
/*  Coulomb stress (cmb) and normal stress (sig)	    */
/************************************************************/

//void cmbfix(double sxx, double syy, double szz, double sxy, double syz, double szx, double p, double f, double *cmb, double *sig, double strike, double dip, double rake) {
//
//
//	double **s, *ns, *ts, *rst, *rdi;
//	double strike_radian, dip_radian, rake_radian, tau, sigg;
//	int i, j;
//
//	s   = dmatrix(1,3,1,3);
//	ns  = dvector(1,3);
//	ts  = dvector(1,3);
//	rst = dvector(1,3);
//	rdi = dvector(1,3);
//
//	s[1][1] = sxx;
//	s[1][2] = sxy;
//	s[1][3] = szx;
//	s[2][1] = sxy;
//	s[2][2] = syy;
//	s[2][3] = syz;
//	s[3][1] = szx;
//	s[3][2] = syz;
//	s[3][3] = szz;
//
//	strike_radian = strike * DEG2RAD;
//	dip_radian    = dip * DEG2RAD;
//	rake_radian   = rake * DEG2RAD;
//
//	ns[1] =  sin(dip_radian) * cos(strike_radian + 0.5 * PI);
//	ns[2] =  sin(dip_radian) * sin(strike_radian + 0.5 * PI);
//	ns[3] = -cos(dip_radian);
//
//	rst[1] = cos(strike_radian);
//	rst[2] = sin(strike_radian);
//	rst[3] = 0.0;
//
//	rdi[1] = cos(dip_radian) * cos(strike_radian + 0.5 * PI);
//	rdi[2] = cos(dip_radian) * sin(strike_radian + 0.5 * PI);
//	rdi[3] = sin(dip_radian);
//
//	/* sigg = normal stress */
//	/*  tau = shear stress  */
//	for (i = 1; i <= 3; i++)
//		ts[i] = rst[i] * cos(rake_radian) - rdi[i] * sin(rake_radian);
//
//	sigg = 0.0;
//	tau = 0.0;
//
//	for (j = 1; j <= 3; j++) {
//		for (i = 1; i <= 3; i++) {
//			sigg = sigg + ns[i] * s[i][j] * ns[j];
//			tau = tau + ts[i] * s[i][j] * ns[j];
//		}
//	}
//
//	/* dCMB = shear stress + friction*(normal stress + pore pressure) */
//	*cmb = tau + f * (sigg + p);
//	*sig = sigg;
//
//	free_dmatrix(s, 1, 3, 1, 3);
//	free_dvector(ns, 1, 3);
//	free_dvector(ts, 1, 3);
//	free_dvector(rst, 1, 3);
//	free_dvector(rdi, 1, 3);
//}



/************************************************/
/*			pscokada						    */
/*												*/
/* x2, y2, z2: receiver coords					*/
/* x1, y1, z1: patch coords						*/
/************************************************/

void pscokada(double x1, double y1, double z1, double strike1, double dip1, double L, double W, double slip_strike, double slip_dip, double open, double x2, double y2, double z2, double *sxx, double *syy, double *szz, double *sxy, double *syz, double *szx, double alpha, double lambda, double mu, double friction) {
//slip, strike in radians.
//x is northward, y is eastward, and z is downward.
//x1, y1, z1 are coordinates of center top of the patch:  -------X-------		(for some reason using center doesn't work).
//														   \             \
//														    \             \
//														     \             \
//															  ---------------

	double depth, stk, sin_stk, cos_stk, sin_2stk, di, csdi, ssdi;
	double DISL1, DISL2, DISL3, AL1, AL2, AW1, AW2, X, Y, Z, UX, UY, UZ, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ;
	double strain1, strain2, strain3, strain4, strain5, strain6, eii, dum1, dum2;
	int IRET;

	// kilometer to meter
	x1 *= KM2M;
	y1 *= KM2M;
	z1 *= KM2M;
	x2 *= KM2M;
	y2 *= KM2M;
	z2 *= KM2M;
	L  *= KM2M;
	W  *= KM2M;

	stk = strike1;
	cos_stk = cos(stk);
	sin_stk = sin(stk);
	sin_2stk = sin(2.0 * stk);

	di = dip1;
	csdi = cos(di);
	ssdi = sin(di);

	/* rake: gegen Uhrzeigersinn  & slip_z: positiv in Richtung Tiefe */
	DISL1 = slip_strike; //cos(DEG2RAD * rake1) * slip;
	DISL2 = slip_dip;    //sin(DEG2RAD * rake1) * slip;
	DISL3 = open;

	// neg und pos Werte? sind das nicht alongdip und againstdip längen?
	AL1 = -0.5 * L;
	AL2 =  0.5 * L;
//	AW1 = -0.5 * W;
//	AW2 =  0.5 * W;
	AW1=-W;		//top of patch is passed to function.
	AW2=0;

	/*  transform from Aki's to Okada's system */
	X = (x2 - x1) * cos_stk + (y2 - y1) * sin_stk;
	Y = (x2 - x1) * sin_stk - (y2 - y1) * cos_stk;
	Z = -z2;
	depth = z1;
	/* z2 corresponds to the recording depth zrec */
	/* z1 corresponds to the depth of slip (reference depth: zref) */

	DC3D(alpha, X, Y, Z, depth, RAD2DEG*di, AL1, AL2, AW1, AW2, DISL1, DISL2, DISL3, &UX, &UY, &UZ, &UXX, &UYX, &UZX, &UXY, &UYY, &UZY, &UXZ, &UYZ, &UZZ, &IRET);

	/* transform from Okada's to Aki's system */
	strain1 = UXX * cos_stk * cos_stk + UYY * sin_stk * sin_stk + 0.5 * (UXY + UYX) * sin_2stk;
	strain2 = UXX * sin_stk * sin_stk + UYY * cos_stk * cos_stk - 0.5 * (UXY + UYX) * sin_2stk;
	strain3 = UZZ;
	strain4 = 0.5 * (UXX - UYY) * sin_2stk - 0.5 * (UXY + UYX) * cos(2.0 * stk);
	dum1 = -0.5 * (UZX + UXZ);
	dum2 = 0.5 * (UYZ + UZY);
	strain5 = dum1 * sin_stk + dum2 * cos_stk;
	strain6 = dum1 * cos_stk - dum2 * sin_stk;
	eii = strain1 + strain2 + strain3;

	dum1 = lambda * eii;
	dum2 = 2.0 * mu;
	*sxx = dum1 + dum2 * strain1;
	*syy = dum1 + dum2 * strain2;
	*szz = dum1 + dum2 * strain3;
	*sxy =        dum2 * strain4;
	*syz =        dum2 * strain5;
	*szx =        dum2 * strain6;

	//p = 0.0;
	//cmbfix(sxx, syy, szz, sxy, syz, szx, p, friction, &cmb, &sig, strike2, dip2, rake2);

	//*DCFS = cmb;
}
