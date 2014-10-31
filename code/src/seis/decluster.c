/*
 * decluster.c
 *
 *  Created on: Nov 11, 2013
 *      Author: camcat
 */

#include "decluster.h"

//int *decluster_catalog_rescalegrid(struct catalog cat, struct crust crst, double Mmain, double **time_missing, int d3){
///* For each grid point, calculates ratio of (tot. time - time revmoved)/(total time), which ca be re used to rescale the rate;
// * but it's better to weight 	remaining earthquakes instead (due to smoothing). Function: decluster_catalog does this.
// * Returns array containing flag (1/0) for selected/excluded events.
// * time_missing is a pointer to a 1d array which will contain the no. of days missing from each grid point (i.i. no. of days inside time window).
// * If time_missing ==NULL, ignored; if *time_missing==NULL, memory will be allocated.
// * d3: flag indicating if 3d distance (instead of horizontal distance) should be used.
// * events must be sorted chronologically.
// */
//	double D, T, d, *tnow, dt;
//	int j;
//	int *sel=ivector(1,cat.Z);
//	for (int i=1; i<=cat.Z; i++) sel[i]=1;
//	int NP= (d3)? crst.N_allP : crst.nLat*crst.nLon;
//
//	if (time_missing){
//		tnow=dvector(1,NP);
//		if (!(*time_missing)) *time_missing=dvector(1,NP);
//		for (int p=1; p<=NP; p++) {
//			tnow[p]=cat.tstart;
//			(*time_missing)[p]= 0.0;
//		}
//	}
//
//	for (int i=1; i<=cat.Z; i++){
//		if (cat.mag[i]>=Mmain){
//			KG74(cat.mag[i], &D, &T);
//			j=i;
//			//decluster catalog:
//			while(j<=cat.Z && (cat.t[j]-cat.t[i])<=T){
//				d= pow(cat.x0[j]-cat.x0[i],2)+pow(cat.y0[j]-cat.y0[i],2);
//				if (d3) d+=pow(cat.depths0[j]-cat.depths0[i],2);
//				d=sqrt(d);
//				if (d<=D) sel[j]=0;
//				j++;
//			}
//			//calculate period missing for each grid point:
//			if (time_missing){
//				for (int p=1; p<=NP; p++){
//					d= pow(crst.x[p]-cat.x0[i],2)+pow(crst.y[p]-cat.y0[i],2);
//					if (d3) d+=pow(crst.depth[p]-cat.depths0[i],2);
//					d=sqrt(d);
//					if (d<=D) {
//						dt=fmin(cat.tend, cat.t[i]+T)-fmax(cat.t[i], tnow[p]);
//						if (dt>0) (*time_missing)[p]+= dt;
//						tnow[p]=fmax(fmin(cat.tend, cat.t[i]+T),tnow[p]);
//					}
//				}
//			}
//		}
//	}
//
//	if (time_missing) free_dvector(tnow, 1, NP);
//	return sel;
//}


int *decluster_catalog(struct catalog cat, double Mmain, double **weights, int d3){
/* Returns array containing flag (1/0) for selected/excluded events.
 * time_missing is a pointer to a 1d array whichweights==NULL, memory will be allocated.
 * d3: flag indicating if 3d distance (instead of horizontal distance) should be used.
 * events must be sorted chronologically.
 */
	double D, T, d, *tnow, *time_missing, dt;
	int *sel=ivector(1,cat.Z);
	for (int i=1; i<=cat.Z; i++) sel[i]=1;

	if (weights){
		tnow=dvector(1,cat.Z);
		time_missing=dvector(1,cat.Z);
		if (!(*weights)) *weights=dvector(1,cat.Z);
		for (int p=1; p<=cat.Z; p++) {
			tnow[p]=cat.tstart;
		}
	}

	for (int i=1; i<=cat.Z; i++) {
		sel[i]=1;
		if (weights) time_missing[i]=0.0;
	}

	for (int i=1; i<=cat.Z; i++){
		if (cat.mag[i]>=Mmain){
			KG74(cat.mag[i], &D, &T);
			for (int j=1; j<=cat.Z; j++){
				if (sel[j]==0) continue;
				d= pow(cat.x0[j]-cat.x0[i],2)+pow(cat.y0[j]-cat.y0[i],2);
				if (d3) d+=pow(cat.depths0[j]-cat.depths0[i],2);
				d=sqrt(d);
				if (d<=D){	//decluster catalog:
					if((cat.t[j]-cat.t[i])>0 && (cat.t[j]-cat.t[i])<=T && j!=i) sel[j]=0;
					else {	//calculate period missing for each grid point:
						dt=fmin(cat.tend, cat.t[i]+T)-fmax(cat.t[i], tnow[j]);
						if (dt>0) time_missing[j]+= dt;
						tnow[j]=fmax(fmin(cat.tend, cat.t[i]+T),tnow[j]);
					}
				}
			}
		}
	}

	if (weights) {
		for (int i=1; i<=cat.Z; i++) (*weights)[i]= (sel[i]==0)? 0.0 : fmin((cat.tend-cat.tstart)/(cat.tend-cat.tstart-time_missing[i]),1.0);
		free_dvector(time_missing, 1, cat.Z);
		free_dvector(tnow, 1, cat.Z);
	}
	return sel;
}

void KG74(double M, double *D, double *T){
/* calculation of Knopoff-Gardner 74 criterium for aftershock declustering
 * function modified from one written by Olga Zakharova.
 *
 * M= magnitude.
 * D= spatial window.
 * T= temporal window.
 */

  if (D) *D = pow(10,(0.1238*M+0.983));
  if (T) *T = (M >= 6.5) ? pow(10,(0.032*M+2.7389)) : pow(10,(0.5409*M-0.547));

  return;
}
