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


void init_crst(struct crust *crst){
/* Initialize variables in crust structure to default values. */


	(*crst).S=NULL;
	(*crst).list_allP=NULL;
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
/* Initialize variables in catalog structure to default values. */

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
	(*cat).weights=dmatrix_firstlevel(Zsel);
	(*cat).b=1.0;

}

struct eqkfm *eqkfm_array(long n1, long n2){
/* Allocate memory to array of eqkfm. */
	struct eqkfm *v;
	v= (struct eqkfm *) malloc((size_t) ((n2-n1+1+NR_END)*sizeof(struct eqkfm)));

	for (int i=NR_END; i<=n2-n1+NR_END; i++){
		v[i].slip_str= NULL;
		v[i].slip_dip= NULL;
                v[i].ts=NULL;
		v[i].nosnap=0;
		v[i].allslip_str= NULL;
                v[i].allslip_dip= NULL;
		v[i].pos_s= NULL;
		v[i].pos_d= NULL;
		v[i].selpoints= NULL;
		v[i].distance= NULL;
		v[i].is_slipmodel=0;
		v[i].np_st=v[i].np_di=v[i].whichfm=v[i].nsel=v[i].noise=0;
		v[i].t=0;
		v[i].lat=0;
		v[i].lon=0;
		v[i].depth=0;
		v[i].mag=0;
		v[i].tot_slip=dvector(0,0);	//only need one element for earthquake sources, will reallocate for afterslip.
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
/* Allocate memory to array of eqkfm. */
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
	return v-n1+NR_END;
}

struct pscmp *pscmp_arrayinit(struct crust v0, long n1, long n2){
/* Allocate memory to array of eqkfm, and initialize variables copied from v0. */

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
		if (v[f].tot_slip) free_dvector(v[f].tot_slip,0,0);
		if (v[f].distance) free_dvector(v[f].distance,1,0);
		if (v[f].pos_s) free_dvector(v[f].pos_s,1,0);
		if (v[f].pos_d) free_dvector(v[f].pos_d,1,0);
		if (v[f].slip_str) free_dvector(v[f].slip_str,1,0);
		if (v[f].slip_dip) free_dvector(v[f].slip_dip,1,0);
		if (v[f].allslip_str) free_dmatrix(v[f].allslip_dip,0,0,1,0);
		if (v[f].allslip_dip) free_dmatrix(v[f].allslip_str,0,0,1,0);

	}
}

// todo should use one of these functions and delete the other one.
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
/* Deallocates memory from variables in catalog structure.
 */

	//assumes that elements have been initialized at position 1.
	free_dvector(cat.t,1, 0);
	free_dvector(cat.mag,1, 0);
	free_dvector(cat.lat0,1, 0);
	free_dvector(cat.lon0,1, 0);
	free_dvector(cat.x0,1, 0);
	free_dvector(cat.y0,1, 0);
	free_dvector(cat.depths0,1, 0);
	free_ivector(cat.ngrid,1, 0);
	free_imatrix_firstlevel(cat.ngridpoints,1,cat.Z,1,0); //uses 0 for upper index, since it doesn't matter (check nrutil.c).
	free_dmatrix_firstlevel(cat.weights,1,cat.Z, 1,0);
}

