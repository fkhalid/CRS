/*
 * mem_mgmt.c
 *
 *Allocates/deallocates memory for various structures.
 *
 *  Created on: Jun 5, 2013
 *      Author: camcat
 */


#include "mem_mgmt.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "../defines.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"

//void shift_cat(struct catalog *cat, int N){
//	/*shifts all vectors in struct cat so than element N becomes the first one.
//	 * NB: shift back (with N->-N+2) before deallocating memory to avoid seg fault.
//	 */
//
//	(*cat).t+=N-1;
//	(*cat).mag+=N-1;
//	(*cat).lat0+=N-1;
//	(*cat).lon0+=N-1;
//	(*cat).x0+=N-1;
//	(*cat).y0+=N-1;
//	(*cat).depths0+=N-1;
//	(*cat).err+=N-1;
//	(*cat).verr+=N-1;
//	(*cat).ngrid+=N-1;
//	(*cat).ngridpoints+=N-1;
//	(*cat).weights+=N-1;
//
//	(*cat).Z-=N-1;
//
//	return;
//}


void init_crst(struct crust *crst){

	(*crst).S=NULL;
	(*crst).list_allP=NULL;		// list. of points (should be same as DCFS0): [1,2,3,...N_allP].
	(*crst).lat=NULL;
	(*crst).lon=NULL;
	(*crst).depth=NULL;
	(*crst).lat_out=NULL;
	(*crst).lon_out=NULL;
	(*crst).depth_out=NULL;
	(*crst).dAgrid=NULL;
	(*crst).str0=dvector(0,0);
	(*crst).dip0=dvector(0,0);
	(*crst).rake0=dvector(0,0);
	(*crst).x=NULL;
	(*crst).y=NULL;
	(*crst).rate0=NULL;
	(*crst).mags=NULL;
	(*crst).GRmags=NULL;
	(*crst).nofmzones=1;
	(*crst).fmzone=NULL;
	(*crst).variable_fixmec=0;

	return;
}

void init_cat1(struct catalog *cat, int Zsel){

	(*cat).Z=Zsel;
	(*cat).t = dvector(1, Zsel);
	(*cat).mag = dvector(1, Zsel);
	(*cat).lat0 = dvector(1, Zsel);
	(*cat).lon0 = dvector(1, Zsel);
	(*cat).x0 = dvector(1, Zsel);
	(*cat).y0 = dvector(1, Zsel);
	(*cat).depths0 = dvector(1, Zsel);
	(*cat).err = dvector(1, Zsel);
	(*cat).verr = dvector(1, Zsel);
	(*cat).ngrid = ivector(1, Zsel);
	//just allocate first level since subarrays may have different length (and will be initialized later).
	(*cat).ngridpoints=imatrix_firstlevel(Zsel);
	(*cat).weights=dmatrix_firstlevel(Zsel);		// weights[0] indicates the fraction of the Gaussian ellipsoid outside the grid.
	(*cat).b=1.0;

}

//void init_cat2(struct catalog *cat, int N, struct crust crst){
//	(*cat).depthmin=crst.depth[1]-0.5*crst.ddepth;
//	(*cat).depthmax=crst.depth[N]+0.5*crst.ddepth;
//	(*cat).latmin=crst.lat[1]-0.5*crst.dlat;
//	(*cat).lonmin=crst.lon[1]-0.5*crst.dlon;
//	(*cat).latmax=crst.lat[N]+0.5*crst.dlat;
//	(*cat).lonmax=crst.lon[N]+0.5*crst.dlon;
//	(*cat).dAeq=pow(Re*PI/180,2)*crst.dlon*crst.dlat;
//	(*cat).latgrid=crst.lat;
//	(*cat).longrid=crst.lon;
//	(*cat).layers=crst.depth;
//}

// todo [coverage] this block is never tested
struct set_of_models *set_of_models_array(long n1, long n2){
/* allocate memory to array of eqkfm. */
	struct set_of_models *v;
	v= (struct set_of_models *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(struct set_of_models)));

	for (int i=NR_END; i<=n2-n1+NR_END; i++){
		v[i].Nmod=1;
		v[i].NF_models=NULL;	//no. of faults for each model;
		v[i].set_of_eqkfm=NULL; //contains all models (size: sum(NF_model)).
	}
	return v-n1+NR_END;
}

struct eqkfm *eqkfm_array(long n1, long n2){
/* allocate memory to array of eqkfm. */
	struct eqkfm *v;
	v= (struct eqkfm *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(struct eqkfm)));

	for (int i=NR_END; i<=n2-n1+NR_END; i++){
		v[i].slip_str= NULL;
		v[i].slip_dip= NULL;
		v[i].pos_s= NULL;
		v[i].pos_d= NULL;
		v[i].selpoints= NULL;
		v[i].distance= NULL;
		v[i].is_slipmodel=0;
		v[i].is_mainshock=0;
		v[i].np_st=v[i].np_di=v[i].whichfm=v[i].nsel=v[i].noise=0;
		v[i].t=0;
		v[i].lat=0;
		v[i].lon=0;
		v[i].depth=0;
		v[i].mag=0;
		v[i].tot_slip=0;
		v[i].L=0;
		v[i].W=0;
		v[i].str1=0;
		v[i].str2=0;
		v[i].dip1=0;
		v[i].dip2=0;
		v[i].rake1=0;
		v[i].rake2=0;
		v[i].index_cat=0;
		v[i].cuts_surf=0;

	}
	return v-n1+NR_END;
}

struct pscmp *pscmp_array(long n1, long n2){
/* allocate memory to array of eqkfm. */
	struct pscmp *v;
	v= (struct pscmp *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(struct pscmp)));

	for (int i=NR_END; i<=n2-n1+NR_END; i++){
		v[i].t=0.0;
//		v[i].dlat=0;
//		v[i].dlon=0;
//		v[i].ddepth=0;
//		v[i].lat=NULL;
//		v[i].lon=NULL;
//		v[i].depth=NULL;
		v[i].fdist=NULL;	//distance to fault
		v[i].S=NULL;
		v[i].cmb=NULL;
		v[i].Z=0;
//		v[i].east_min=0;
//		v[i].east_max=0;
//		v[i].north_min=0;
//		v[i].north_max=0;
		v[i].st1=NULL;
		v[i].di1=NULL;
		v[i].ra1=NULL;
		v[i].st2=NULL;
		v[i].di2=NULL;
		v[i].ra2=NULL;
		v[i].which_pts=NULL;
		//v[i].nLat=v[i].nLon=v[i].nD=0;	//describing overall geometry.
		v[i].nsel=0;
	    v[i].index_cat=0;
	    v[i].NF=0;
	}
	return v-n1+NR_END;
}

struct pscmp *pscmp_arrayinit(struct crust v0, long n1, long n2){
/* allocate memory to array of eqkfm, and initialize variables copied from v0. */
		struct pscmp *v;
		v= (struct pscmp *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(struct pscmp)));

		for (int i=NR_END; i<=n2-n1+NR_END; i++){
			v[i].t=0.0;
			v[i].fdist=NULL;	//distance to fault
			v[i].S=NULL;
			v[i].cmb=NULL;
			v[i].Z=0;
			v[i].st1=NULL;
			v[i].di1=NULL;
			v[i].ra1=NULL;
			v[i].st2=NULL;
			v[i].di2=NULL;
			v[i].ra2=NULL;
			v[i].which_pts=NULL;
			v[i].nsel=0;
		    v[i].index_cat=0;
		    v[i].NF=0;
		}
		return v-n1+NR_END;;
}

// todo [coverage] this block is never tested
void free_eqkfmarray(struct eqkfm *v, long n1, long n2){
/* free a eqkfm vector allocated with eqkfm_array() */
	free((FREE_ARG) (v+n1-NR_END));
}

// todo [coverage] this block is never tested
void freefull_eqkfmarray(struct eqkfm *v, long n1, long n2){

	for (int f=n1; f<=n2; f++){
		if (v[f].distance) free(v[f].distance);
		if (v[f].pos_s) free(v[f].pos_s);
		if (v[f].pos_d) free(v[f].pos_d);
		if (v[f].slip_str) free(v[f].slip_str);
		if (v[f].slip_dip) free(v[f].slip_dip);
	}
	//free((FREE_ARG) (v+n1-NR_END));
}

// todo [coverage] this block is never tested
void freepart_pscmparray(struct pscmp *v, long n1, long n2){
//only frees stuff which wasn't linked to other structure (see pscmp_arrayinit).
	for (int i=n1; i<=n2; i++){

		if (v[i].S!=NULL) {
			for (int j=1; j<=v[i].nsel; j++) free_dmatrix(v[i].S[j],1,3,1,3);
			free(v[i].S[i]);
		}
		if (v[i].cmb!=NULL) free(v[i].cmb);
		if (v[i].st1!=NULL) free(v[i].st1);
		if (v[i].di1!=NULL) free(v[i].di1);
		if (v[i].ra1!=NULL) free(v[i].ra1);
		if (v[i].st2!=NULL) free(v[i].st2);
		if (v[i].di2!=NULL) free(v[i].di2);
		if (v[i].ra2!=NULL) free(v[i].ra2);
	}
	free((FREE_ARG) (v+n1-NR_END));
}

// todo [coverage] this block is never tested
void freefull_pscmparray(struct pscmp *v, long n1, long n2){

	for (int i=n1; i<=n2; i++){

//		if (v[i].lat!=NULL) free(v[i].lat);
//		if (v[i].lon!=NULL) free(v[i].lon);
//		if (v[i].depth!=NULL) free(v[i].depth);
		if (v[i].fdist!=NULL) free(v[i].fdist);
		if (v[i].S!=NULL) {
			for (int j=1; j<=v[i].nsel; j++) free_dmatrix(v[i].S[j],1,3,1,3);
			free(v[i].S[i]);
		}
		if (v[i].cmb!=NULL) free(v[i].cmb);
		if (v[i].st1!=NULL) free(v[i].st1);
		if (v[i].di1!=NULL) free(v[i].di1);
		if (v[i].ra1!=NULL) free(v[i].ra1);
		if (v[i].st2!=NULL) free(v[i].st2);
		if (v[i].di2!=NULL) free(v[i].di2);
		if (v[i].ra2!=NULL) free(v[i].ra2);
		if (v[i].which_pts!=NULL) free(v[i].which_pts);
	}
	free((FREE_ARG) (v+n1-NR_END));
}

void free_cat(struct catalog cat){
//assumes that elements have been initialized at position 1.
//uses 0 for upper index, since it doesn't matter (check nrutil.c).

	free_dvector(cat.t,1, 0);
	free_dvector(cat.mag,1, 0);
	free_dvector(cat.lat0,1, 0);
	free_dvector(cat.lon0,1, 0);
	free_dvector(cat.x0,1, 0);
	free_dvector(cat.y0,1, 0);
	free_dvector(cat.depths0,1, 0);
	free_ivector(cat.ngrid,1, 0);
	free_imatrix_firstlevel(cat.ngridpoints,1,cat.Z,1,0);
	free_dmatrix_firstlevel(cat.weights,1,cat.Z, 1,0);
//	free_dvector(cat.xgrid,1,0);
//	free_dvector(cat.ygrid,1,0);
//	free_dvector(cat.dAgrid,1,0);
}

