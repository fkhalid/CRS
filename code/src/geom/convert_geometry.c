/*
 * convert_geometry.c
 *
 *  Created on: Sep 13, 2013
 *      Author: camcat
 */

#include "convert_geometry.h"

int convert_geometry(struct crust crst, double *old_v, double **new_v, int sum, int increase_resolution){
// if increase_resolution, will use crst.nX_out as starting points, and calculate new values for crst.nX; otherwise, will convert from larger to smaller cells.
// for sum==0, quantity is an intensive variable -> take average. for sum==1, takes sum (only for increase_resolution==0).
// points change along lat, then along lon, then along depth.
// indices start from 1.
// *new_v will be initialized if NULL; otherwise, must have correct no. of elements!

	int P[3], Pn[3], nsub[3];
	int D1, D2, D3, D1D2;	//no o fpoints per dimension (old geometry);
	int D1n, D2n, D3n;	//no o fpoints per dimension (new geometry);
	int NP, NPn, nsub_tot;
	int ind;

	NPn= (increase_resolution)? crst.N_allP : crst.nLat_out*crst.nLon_out*crst.nD_out;
	NP= (increase_resolution)? crst.nLat_out*crst.nLon_out*crst.nD_out : crst.N_allP;

	if (!crst.uniform || NP==NPn){
		*new_v=old_v;
		return 0;
	}

	if (increase_resolution){
		D1=crst.nLat_out;
		D2=crst.nLon_out;
		D3=crst.nD_out;
		D1D2=D1*D2;
		D1n=crst.nLat;
		D2n=crst.nLon;
		D3n=crst.nD;
	}

	else {
		D1=crst.nLat;
		D2=crst.nLon;
		D3=crst.nD;
		D1D2=D1*D2;
		D1n=crst.nLat_out;
		D2n=crst.nLon_out;
		D3n=crst.nD_out;
	}

	//no. of grid points between high resolution geometry and low resolution geometry (per dimension);
	nsub[0]=crst.nLat/crst.nLat_out;
	nsub[1]=crst.nLon/crst.nLon_out;
	nsub[2]=crst.nD/crst.nD_out;
	nsub_tot=nsub[0]*nsub[1]*nsub[2];

	if (crst.nLat%crst.nLat_out!=0 || crst.nLon%crst.nLon_out!=0 || crst.nD%crst.nD_out!=0) {
		if (verbose_level>0) printf(" ** Error: calculation cells are not a multiple of forecast cells - can not recalculate geometry -> Using old geometry. (convert_geometry)\n");
		*new_v=old_v;
		return(1);
	}

	if (!(*new_v)) *new_v=dvector(1,NPn);
	for (int i=1; i<=NPn; i++) (*new_v)[i]=0.0;

	for (int i=1; i<=NP; i++){
		//reshape linear array into 3x3 array: i -> (P1,P2,P3).
		P[0]=(i-1)%D1+1;
		P[1]=((i-1)%D1D2)/D1 +1;
		P[2]=(i-1)/D1D2 +1;

		if (increase_resolution){
			for (int x=1; x<=nsub[0]; x++){
				Pn[0]=(P[0]-1)*nsub[0]+x;
				for (int y=1; y<=nsub[1]; y++){
					Pn[1]=(P[1]-1)*nsub[1]+y;
					for (int z=1; z<=nsub[2]; z++){
						Pn[2]=(P[2]-1)*nsub[2]+z;
						ind=Pn[0]+D1n*(Pn[1]-1)+D1n*D2n*(Pn[2]-1);
						(*new_v)[ind]= old_v[i];
					}
				}
			}
		}

		else {
			//calculate new indices:
			for (int n=0; n<3; n++) Pn[n]=(P[n]-1)/nsub[n]+1;

			//calculate linear index:
			ind=Pn[0]+D1n*(Pn[1]-1)+D1n*D2n*(Pn[2]-1);
			(*new_v)[ind]= (sum)? (*new_v)[ind]+old_v[i] : (*new_v)[ind]+old_v[i]*(1.0/nsub_tot);
		}
	}

	return(0);

}




