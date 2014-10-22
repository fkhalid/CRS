///*
// * dist2fault.c
// *
// *  Created on: Apr 18, 2013
// *      Author: camcat
// */
//
//#include "dist2fault.h"
//
//
//void dist2faultDCFS(struct pscmp *DCFS, struct crust crst, struct eqkfm *eqkfm1){
//
//	(*DCFS).fdist=dvector(1,(*DCFS).nsel);	//TODO deallocate.
//	dist2fault(crst.lat, crst.lon, crst.depth, (*DCFS).fdist, (*DCFS).nsel, (*DCFS).NF, eqkfm1);
//
//}
//
//double *dist2fault0(double *lats, double *lons, double *depths, int NP, double strike0, double dip0,
//		double Lat0, double Lon0, double D0, double *pos_s, double *pos_d){
///* Translated from matlab function with same name.
// * assumes rectangular fault.
// *
// * input:
// * lats, lons, depths: vectors [1...NP], containing coordinates of grid points.
// * Lat0, Lon0, D0: coordinates of a reference point on the fault plane.
// * pos_s[0,1]: the distance along strike from the reference point to the sides edges of the fault.
// * pos_d[0,1]: the distance along dip from the reference point to the top/bottom edged of the fault.
// *
// * output:
// * vector with indices [1...NP] containing distance (positive or negative depending on opposite sides of the fault).
// */
//
//double *dist;
//double *dipV, *A, *A0, *B, *B0, *P, *J, *D;
//double l, d, dir, dz1, dz2, Dmin, Dmax, Depth;
//double strike, dip;
//
//strike=strike0*DEG2RAD;
//dip=dip0*DEG2RAD;
//
//dist=dvector(1,NP);
//dipV=dvector(1,3);
//A=dvector(1,3);
//A0=dvector(1,3);
//B=dvector(1,3);
//B0=dvector(1,3);
//J=dvector(1,2);
//P=dvector(1,3);
//D=dvector(1,3);
//
//	// top/bottom depth.
//	dz1=sin(dip)*pos_d[0];
//	dz2=sin(dip)*pos_d[1];
//	Dmin=D0+dz1;
//	Dmax=D0+dz2;
//
//	// coordinates of top corners:
//	A0[1]=pos_s[0]*sin(strike)+pos_d[0]*cos(dip)*cos(strike);	//x coordinate: east.
//	A0[2]=pos_s[0]*cos(strike)-pos_d[0]*cos(dip)*sin(strike);
//	A0[3]=Dmin;
//	B0[1]=pos_s[1]*sin(strike)+pos_d[0]*cos(dip)*cos(strike);
//	B0[2]=pos_s[1]*cos(strike)-pos_d[0]*cos(dip)*sin(strike);
//	B0[3]=Dmin;
//
//	//dip vector
//	dipV[1]=cos(dip)*sin(strike+PI/2.0);
//	dipV[2]=cos(dip)*cos(strike+PI/2.0);
//	dipV[3]=sin(dip);
//
//	for (int k=1; k<=NP; k++){
//		P[1]=Re*(lons[k]-Lon0)*PI/180*cos(Lat0*PI/180);
//		P[2]=Re*(lats[k]-Lat0)*PI/180;
//		P[3]=depths[k];
//
//		dist2line(P,A0,B0,2,0,J, &l);
//		dir=vdotv(J,dipV,2);
//		l=-l*(dir/fabs(dir));
//
//		Depth=P[3]-Dmin+l*sin(dip)*cos(dip)-(P[3]-Dmin)*cos(dip)*cos(dip);
//
//		if (Depth>Dmax-Dmin) Depth=Dmax-Dmin;
//		else if (Depth<0) Depth=0;
//
//		for (int j=1; j<=3; j++) {		//A, B become coordinates of points at correct depth.
//			A[j]= A0[j]+(Depth/sin(dip))*dipV[j];
//			B[j]= B0[j]+(Depth/sin(dip))*dipV[j];
//		}
//
//		dist2line(P,A,B,3,1,D, &d);
//
//		dir=vdotv(D,dipV,2);
//		dist[k]=-d*(dir/fabs(dir));	//wedge: +ve, below: -ve. (may be different in Matlab version).
//
//	}
//
//	return dist;
//}
//
//
//
//void dist2fault(double *lats, double *lons, double *depths, double *dist, int NP, int NF, struct eqkfm *eqkfm1){
//// Translated from matlab function with same name.
//
//double strike, dip, rake;
//double *dipV, *A, *A0, *B, *B0, *P, *J, *D;
//double l, d, dir, dd, dz1, dz2, Dmin, Dmax, Depth, dpos;
//double Lat0, Lon0, D0;
//
//dipV=dvector(1,3);
//A=dvector(1,3);
//A0=dvector(1,3);
//B=dvector(1,3);
//B0=dvector(1,3);
//J=dvector(1,2);
//P=dvector(1,3);
//D=dvector(1,3);
//
//for (int h=0; h<NF; h++){
//
//	switch (eqkfm1[h].whichfm){
//		case 1:
//			strike=DEG2RAD*eqkfm1[h].str1;
//			dip=DEG2RAD*eqkfm1[h].dip1;
//			rake=DEG2RAD*eqkfm1[h].rake1;
//			break;
//		case 2:
//			strike=DEG2RAD*eqkfm1[h].str2;
//			dip=DEG2RAD*eqkfm1[h].dip2;
//			rake=DEG2RAD*eqkfm1[h].rake2;
//	}
//
//	Lat0=eqkfm1[h].lat;
//	Lon0=eqkfm1[h].lon;
//	D0=eqkfm1[h].depth;
//
//	// top/bottom depth.
//	dd=eqkfm1[h].W/eqkfm1[h].np_di;
//	dz1=sin(dip)*eqkfm1[h].pos_d[1];
//	dz2=sin(dip)*(eqkfm1[h].pos_d[eqkfm1[h].np_di*eqkfm1[h].np_st]+dd);
//	Dmin=D0+dz1;
//	Dmax=D0+dz2;
//
//	// coordinates of top corners:
//	int fp=1, lp=1;
//	for (int c=1; c<=eqkfm1[h].np_di*eqkfm1[h].np_st; c++){
//		if (eqkfm1[h].pos_d[c]==eqkfm1[h].pos_d[1]){
//			if (eqkfm1[h].pos_s[c]<eqkfm1[h].pos_s[fp]) fp=c;
//			if (eqkfm1[h].pos_s[c]>eqkfm1[h].pos_s[lp]) lp=c;
//		}
//	}
//	dpos=eqkfm1[h].L/eqkfm1[h].np_st;
//	A0[1]=(eqkfm1[h].pos_s[fp]-0.5*dpos)*sin(strike)+eqkfm1[h].pos_d[fp]*cos(dip)*cos(strike);	//x coordinate: east.
//	A0[2]=(eqkfm1[h].pos_s[fp]-0.5*dpos)*cos(strike)-eqkfm1[h].pos_d[fp]*cos(dip)*sin(strike);
//	A0[3]=Dmin;
//	B0[1]=(eqkfm1[h].pos_s[lp]+0.5*dpos)*sin(strike)+eqkfm1[h].pos_d[lp]*cos(dip)*cos(strike);
//	B0[2]=(eqkfm1[h].pos_s[lp]+0.5*dpos)*cos(strike)-eqkfm1[h].pos_d[lp]*cos(dip)*sin(strike);
//	B0[3]=Dmin;
//
//	//dip vector
//	dipV[1]=cos(dip)*sin(strike+PI/2.0);
//	dipV[2]=cos(dip)*cos(strike+PI/2.0);
//	dipV[3]=sin(dip);
//
//	//dlat=(180*dy/pi)/R;
//	//dlon=(180*dx/pi)/(R*cos(Lat0*pi/180));
//	//Lats=Lat0+dlat;
//	//Lons=Lon0+dlon;
//	//Ds=D0+dz;
//	//Dmin=D0+min(dz(:));
//	//Dmax=D0+max(dz(:));
//
//	for (int k=1; k<=NP; k++){
//
//		//latlon2localcartesian((*DCFS).lat[k],(*DCFS).lon[k], Lat0, Lon0, P+2, P+1);
//		P[1]=Re*(lons[k]-Lon0)*PI/180*cos(Lat0*PI/180);
//		P[2]=Re*(lats[k]-Lat0)*PI/180;
//		P[3]=depths[k];
//
//		dist2line(P,A0,B0,2,0,J, &l);	//determine direction.
//		dir=vdotv(J,dipV,2);
//		l=-l*(dir/fabs(dir));
//
//		Depth=P[3]-Dmin+l*sin(dip)*cos(dip)-(P[3]-Dmin)*cos(dip)*cos(dip);
//
//		if (Depth>Dmax-Dmin) Depth=Dmax-Dmin;
//		else if (Depth<0) Depth=0;
//
//		for (int j=1; j<=3; j++) {		//A, B become coordinates of points at correct depth.
//			A[j]= A0[j]+(Depth/sin(dip))*dipV[j];
//			B[j]= B0[j]+(Depth/sin(dip))*dipV[j];
//		}
//
//		dist2line(P,A,B,3,1,D, &d);	//TODO maybe set finiteline=1 (currently assumes fault continues forever).
//
//		if (h==0 || fabs(d)<fabs(dist[k])){
//			dir=vdotv(D,dipV,2);
//			dist[k]=-d*(dir/fabs(dir));	//wedge: +ve, below: -ve. (may be different in Matlab version).
//			//Vdist(k,:)=D;
//		}
//	}
//}
//}
//
//void dist2line(double *P, double *A, double *B, int d3, int finitesegment, double *D, double *d){
//// Translated from matlab function with same name.
//
//	double a,b,l, ap, bp, cosalpha, cal;
//	double *L, *AP, *BP;
//
//	L=dvector(1,d3);
//	AP=dvector(1,d3);
//	BP=dvector(1,d3);
//
//	a=norm(A,d3);
//	b=norm(B,d3);
//	for (int k=1; k<=d3; k++) {
//		L[k]=B[k]-A[k];
//		AP[k]=A[k]-P[k];
//		BP[k]=B[k]-P[k];
//	}
//
//	l=norm(L,d3);
//	ap=norm(AP,d3);
//	bp=norm(BP,d3);
//
//	cal=vdotv(AP,L,d3);
//	cosalpha=cal/(ap*l);
//
//	for (int k=1; k<=d3; k++) D[k]=A[k]-L[k]*(ap*cosalpha/l) - P[k];
//	*d=norm(D,d3);
//
//	if (finitesegment==1){
//		if ((pow(ap,2)-pow(*d,2))>pow(l,2) || (pow(bp,2)-pow(*d,2))>pow(l,2)){
//			for (int k=1; k<=d3; k++) D[k]= (ap<bp)? A[k]-P[k] : B[k]-P[k];
//			*d=norm(D,d3);
//		}
//	}
//}
//
//void select_onofffault(double *FMdist, int NFM, int **labove, int **lbelow, int **lon, int *Nab, int *Nbe, int *Non){
//
//	int nabove,nbelow,non;
//	nabove=nbelow=non=0;
//	*Nab=*Nbe=*Non=0;
//	for (int k=1; k<=NFM; k++){
//		if (FMdist[k]>=30) (*Nab)++;
//		else {
//			if (FMdist[k]<=-30) (*Nbe)++;
//			else (*Non)++;
//		}
//	}
//
//	(*labove)=ivector(1,*Nab);	//TODO deallocate at the end.
//	(*lbelow)=ivector(1,*Nbe);
//	(*lon)=ivector(1,*Non);
//
//	for (int k=1; k<=NFM; k++){
//		if (FMdist[k]>=30) {
//			nabove+=1;
//			(*labove)[nabove]=k;
//		}
//		else {
//			if (FMdist[k]<=-30) {
//				nbelow+=1;
//				(*lbelow)[nbelow]=k;
//			}
//			else {
//				non+=1;
//				(*lon)[non]=k;
//			}
//		}
//	}
//}
