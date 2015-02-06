/*
 * eqkfm_copy.c
 *
 *  Created on: Nov 20, 2013
 *      Author: camcat
 */

#include "eqkfm_copy.h"

// todo [coverage] this block is never tested
void empty_eqkfm(struct eqkfm *eqkfm0){
	(*eqkfm0).np_st=(*eqkfm0).np_di=0;
	(*eqkfm0).nsel=0;
	(*eqkfm0).nosnap=0;
	(*eqkfm0).tot_slip=dvector(0,0);	//only need one element for earthquake slip models, will reallocate for afterslip.
	(*eqkfm0).L=0;
	(*eqkfm0).W=0;
	(*eqkfm0).ts=NULL;
	(*eqkfm0).tevol=NULL;
	(*eqkfm0).slip_str=NULL;
	(*eqkfm0).slip_dip=NULL;
	(*eqkfm0).allslip_str=NULL;
	(*eqkfm0).allslip_dip=NULL;
	(*eqkfm0).pos_s=NULL;
	(*eqkfm0).pos_d=NULL;
	(*eqkfm0).distance=NULL;
	(*eqkfm0).selpoints=NULL;		//indices of cell points affected by this event.
}

// todo [coverage] this block is never tested
void copy_eqkfm_nolocation_noindex_notime(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
//slightly wasteful, but at least don't need to change it if new variables are added to structure (junk change copy_eqkfm_all).

	double lat, lon, depth, t;
	int i;

	i=(*eqkfm2).index_cat;
	lat=(*eqkfm2).lat;
	lon=(*eqkfm2).lon;
	depth=(*eqkfm2).depth;
	t=(*eqkfm2).t;

	copy_eqkfm_all(eqkfm1, eqkfm2);

	(*eqkfm2).t=t;
	(*eqkfm2).index_cat=i;
	(*eqkfm2).lat=lat;
	(*eqkfm2).lon=lon;
	(*eqkfm2).depth=depth;

}

void copy_eqkfm_noindex_notime(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
//slightly wasteful, but at least don't need to change it if new variables are added to structure (junk change copy_eqkfm_all).

	double lat, lon, depth, t;
	int i;

	i=(*eqkfm2).index_cat;
	t=(*eqkfm2).t;

	copy_eqkfm_all(eqkfm1, eqkfm2);

	(*eqkfm2).t=t;
	(*eqkfm2).index_cat=i;
}

void copy_eqkfm_attributes(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){

	//general properties:
	(*eqkfm2).is_slipmodel=eqkfm1.is_slipmodel;
    (*eqkfm2).index_cat=eqkfm1.index_cat;	//index of event in catalog (only for the catalog used for LL calculation).
	(*eqkfm2).nsel=eqkfm1.nsel;
	(*eqkfm2).selpoints= eqkfm1.selpoints;
	(*eqkfm2).distance= eqkfm1.distance;
	(*eqkfm2).noise=eqkfm1.noise;
	(*eqkfm2).cuts_surf=eqkfm1.cuts_surf;
	(*eqkfm2).top=eqkfm1.top;
	(*eqkfm2).nosnap=eqkfm1.nosnap;


	//earthquake/afterslip properties:
	(*eqkfm2).t=eqkfm1.t;
	(*eqkfm2).ts=eqkfm1.ts;
	(*eqkfm2).tevol=eqkfm1.tevol;
	(*eqkfm2).mag=eqkfm1.mag;
	(*eqkfm2).lat=eqkfm1.lat;
	(*eqkfm2).lon=eqkfm1.lon;
	(*eqkfm2).x=eqkfm1.x;
	(*eqkfm2).y=eqkfm1.y;
	(*eqkfm2).depth=eqkfm1.depth;

	return;

}

void copy_eqkfm_focmec(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){

	//focal mechanism:
	(*eqkfm2).whichfm=eqkfm1.whichfm;
	(*eqkfm2).str1=eqkfm1.str1;
	(*eqkfm2).dip1=eqkfm1.dip1;
	(*eqkfm2).rake1=eqkfm1.rake1;
	(*eqkfm2).str2=eqkfm1.str2;
	(*eqkfm2).dip2=eqkfm1.dip2;
	(*eqkfm2).rake2=eqkfm1.rake2;

	return;

}

void copy_eqkfm_slipmodel(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){

	//slipmodel properties:
	(*eqkfm2).tot_slip=eqkfm1.tot_slip;
	(*eqkfm2).L=eqkfm1.L;
	(*eqkfm2).W=eqkfm1.W;
	(*eqkfm2).np_st=eqkfm1.np_st;
	(*eqkfm2).np_di=eqkfm1.np_di;
	(*eqkfm2).pos_s=eqkfm1.pos_s;
	(*eqkfm2).pos_d=eqkfm1.pos_d;
	(*eqkfm2).slip_str=eqkfm1.slip_str;
	(*eqkfm2).slip_dip=eqkfm1.slip_dip;
	(*eqkfm2).allslip_str=eqkfm1.allslip_str;
	(*eqkfm2).allslip_dip=eqkfm1.allslip_dip;

	return;

}

void copy_eqkfm_all(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){

	copy_eqkfm_attributes(eqkfm1, eqkfm2);
	copy_eqkfm_focmec(eqkfm1, eqkfm2);
	copy_eqkfm_slipmodel(eqkfm1, eqkfm2);

	return;
}


// todo [coverage] this block is never tested
void copy_eqkfm_noslipmodel(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){

	copy_eqkfm_attributes(eqkfm1, eqkfm2);
	copy_eqkfm_focmec(eqkfm1, eqkfm2);

	return;

}
