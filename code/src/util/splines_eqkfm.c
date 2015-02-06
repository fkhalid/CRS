/*
 * splines_eqkfm.c
 *
 *  Created on: May 16, 2013
 *      Author: camcat
 */

#include <math.h>

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "fit_splines.h"
#include "nrutil.h"

void splines_eqkfm(struct eqkfm **eqkfm_aft, int Nas, int NF, double *times1, double *times2, int L, long *seed){
	double **slipbefore_st, **slipbefore_di, **slip_after_st, **slip_after_di;
	int *NP, NP_tot;

	NP=ivector(0,NF-1);
	NP_tot=0;

	for (int f=0; f<NF; f++) {
		NP[f]= (*eqkfm_aft)[f].np_di*(*eqkfm_aft)[f].np_st;
		NP_tot=fmax(NP_tot,NP[f]);
	}

	slipbefore_st=dmatrix(1,NP_tot,1,Nas);
	slipbefore_di=dmatrix(1,NP_tot,1,Nas);
	slip_after_st=dmatrix(1,NP_tot,1,L);
	slip_after_di=dmatrix(1,NP_tot,1,L);

	for (int f=0; f<NF; f++){
		for (int pt=1; pt<=NP[f]; pt++){
			for (int t=0; t<Nas; t++){
				slipbefore_st[pt][t+1]=(*eqkfm_aft)[f].allslip_str[t][pt];
				slipbefore_di[pt][t+1]=(*eqkfm_aft)[f].allslip_dip[t][pt];
			}
		}
		fit_splines(times1, times2, Nas, L, NP[f], slipbefore_st, (double *) 0, &slip_after_st, 2, seed);
		fit_splines(times1, times2, Nas, L, NP[f], slipbefore_di, (double *) 0, &slip_after_di, 2, seed);

		(*eqkfm_aft)[f].ts=times2;
		(*eqkfm_aft)[f].nosnap=L;
		//reallocate memory:
		free_dmatrix((*eqkfm_aft)[f].allslip_str,0,Nas-1, 1,NP[f]);
		free_dmatrix((*eqkfm_aft)[f].allslip_dip,0,Nas-1, 1,NP[f]);
		free_dvector((*eqkfm_aft)[f].tot_slip,0,Nas-1);
		(*eqkfm_aft)[f].allslip_str=dmatrix(0,L-1,1,NP[f]);
		(*eqkfm_aft)[f].allslip_dip=dmatrix(0,L-1,1,NP[f]);
		(*eqkfm_aft)[f].tot_slip=dvector(0,L-1);

		for (int t=0; t<L; t++){
			for (int pt=1; pt<=NP[f]; pt++){
				(*eqkfm_aft)[f].allslip_str[t][pt]=slip_after_st[pt][t+1];
				(*eqkfm_aft)[f].allslip_dip[t][pt]=slip_after_di[pt][t+1];
			}
		}
	}

	free_ivector(NP,0,NF-1);
	free_dmatrix(slipbefore_st,1,NP_tot,1,Nas);
	free_dmatrix(slipbefore_di,1,NP_tot,1,Nas);
}


//void splines_eqkfm(struct eqkfm *eqkfm_aft, int Nas, int NF, double *times1, double *times2, int L, struct eqkfm *all_afterslips, long *seed){
//	double **slipbefore_st, **slipbefore_di, **slip_after_st, **slip_after_di;
//	int *NP, NP_tot;
//	int pure_strslip, pure_thrustnorm;
//
//	for (int f=1; f<=NF*L; f++) copy_eqkfm_all(eqkfm_aft[1], all_afterslips+f);
//
//	NP=ivector(1,NF);
//	NP_tot=0;
//
//	for (int f=1; f<=NF; f++) {
//		NP[f]= eqkfm_aft[f].np_di*eqkfm_aft[f].np_st;
//		NP_tot=fmax(NP_tot,NP[f]);
//	}
//
//	slipbefore_st=dmatrix(1,NP_tot,1,Nas);
//	slipbefore_di=dmatrix(1,NP_tot,1,Nas);
//	slip_after_st=dmatrix(1,NP_tot,1,L);
//	slip_after_di=dmatrix(1,NP_tot,1,L);
//
//	for (int f=1; f<=NF; f++){
//		pure_thrustnorm=pure_strslip=0;
////		pure_thrustnorm=pure_strslip=1;
////		if (fmod(eqkfm_aft[f].rake1+90,180.0)!=0) pure_thrustnorm=0;
////		if (fmod(eqkfm_aft[f].rake1,180.0)!=0) pure_strslip=0;
//		for (int pt=1; pt<=NP[f]; pt++){
//			for (int t=1; t<=Nas; t++){
//				slipbefore_st[pt][t]=eqkfm_aft[NF*(t-1)+f].slip_str[pt];
//				slipbefore_di[pt][t]=eqkfm_aft[NF*(t-1)+f].slip_dip[pt];
//			}
//		}
//		if (pure_thrustnorm==0) fit_splines(times1, times2, Nas, L, NP[f], slipbefore_st, (double *) 0, &slip_after_st, 2, seed);
//		// todo [coverage] this block is never tested
//		else for (int pt=1; pt<=NP[f]; pt++) for (int t=1; t<=L; t++) slip_after_st[pt][t]=0.0;
//		if (pure_strslip==0) fit_splines(times1, times2, Nas, L, NP[f], slipbefore_di, (double *) 0, &slip_after_di, 2, seed);
//		// todo [coverage] this block is never tested
//		else for (int pt=1; pt<=NP[f]; pt++) for (int t=1; t<=L; t++) slip_after_di[pt][t]=0.0;
//		for (int t=1; t<=L; t++){
//			copy_eqkfm_all(eqkfm_aft[f],all_afterslips+NF*(t-1)+f);
//			all_afterslips[NF*(t-1)+f].slip_str=dvector(1,NP[f]);
//			all_afterslips[NF*(t-1)+f].slip_dip=dvector(1,NP[f]);	//TODO deallocate this, and here also allocate/copy other stuff in structure.
//			for (int pt=1; pt<=NP[f]; pt++){
//				all_afterslips[NF*(t-1)+f].slip_str[pt]=slip_after_st[pt][t];
//				all_afterslips[NF*(t-1)+f].slip_dip[pt]=slip_after_di[pt][t];
//			}
//		}
//	}
//
//	free_ivector(NP,1,NF);
//	free_dmatrix(slipbefore_st,1,NP_tot,1,Nas);
//	free_dmatrix(slipbefore_di,1,NP_tot,1,Nas);
//}
