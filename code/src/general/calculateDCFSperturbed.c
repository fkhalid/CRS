/*
 * calculateDCFSperturbed.c
 *
 *  Created on: Mar 19, 2013
 *      Author: camcat
 */

#include "calculateDCFSperturbed.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "../defines.h"
#include "../okada/okadaDCFS.h"
#include "../seis/cmbopt.h"
#include "../util/error.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "../util/ran1.h"
#include "../inp_out/write_csep_forecast.h"
#include "mem_mgmt.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

void calculateDCFSperturbed(double **DCFSrand, struct pscmp *DCFS, struct eqkfm *eqkfmAf,
							struct eqkfm *eqkfm0, struct flags flag,
							double *times, int Nmain, int NA, struct crust crst,
							struct Coeff_LinkList *AllCoeff, int NTScont,
							double **focmec, int *fmzoneslim, int NFM, long *seed,
							double tdata0, double tdata1, int refresh, int which_recfault) {

/* Input:
 * 
 * eqkfmAf contains afterslip snapshots;
 * eqkfm0 contains seismic sources;
 * times: times to which elements of tevol correspond: S(t=times[j])=S0*tevol[j], where S=slip, S0 is the slip contained in eqkfmAf.
 * NTScont is the total number of time steps for continuous process;
 * Nmain: length of eqkfm0;	//todo change variable name
 * NA: length of eqkfmAf;
 * AllCoeff: okada coefficients for mainshocks (events in eqkfm0);
 * focmec contains sample of focal mechanisms, NFM is its length (1->NFM)
 * fmzoneslin: indices of focmec corresponding to limits of distinct foc. mec. areas;
 * tdata0, tdata1: start and end time for which data should be used;
 * refresh: flag to be set to 1 if the slip models have changed from previous function call;
 * flag contains various flags, including:
 * 		afterslip and aftershocks are flags indicating if these processes should be included.
 * 		vary_recfault is a flag indicating which received faults to use. 0: use fix planes; 1: vary foc. mec. (use which_fm, or if which_fm==0, choose random one); 2: use OOPs.
 * which_recfault gives index of foc.mec to select; if set to 0, choose random one (useful if NFM>>nsur). if vary_recfault=0, this is meaningless.
 *
 *
 * Output:
 *
 * DCFSrand[i][j] contains the ith stress change at gridpoint j due to continuous processes (modeled linearly between time steps).
 * DCFS[k].cmp[j] contains the stress change due to kth event at gridpoint j (modeled as a step).
 *
 *
 */

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	//flags:
	int	afterslip=flag.afterslip, \
		vary_recfault=flag.err_recfault, \
		gridpoints_err=flag.err_gridpoints, \
		multisnap=flag.aseismic_multisnap, \
		full_field=(flag.sources_without_focmec==2);

	static double *strike0, *dip0, *rake0;	//strike0, dip0, rake0 are focal mechanism that will change at each iteration.
	double slip;
	double ***Stemp;
	static double lat0, lon0;
	int rand;
	int NgridT=crst.N_allP;
	int NFsofar=0;
	int last, first, i;
	int a_ev=0;	//counter for DCFS_Af
	int nfaults=0;	//counter for eqkfmAf
	int NTSeff;
	static int fm_offset=0;	//offset to be added to crst.str0 (dip0) if a fixed receiver fault is used (more details below).
	static struct eqkfm *eqkfm2;
	static struct eqkfm *eqkfm2A;
	struct Coeff_LinkList *temp, *AllCoeffaft;
	float ***Coeffs_st, ***Coeffs_dip, ***Coeffs_open;	//coefficients for tensor;
	static float **Coeff_ResS, **Coeff_ResD;	//coefficients for cmb (resolved).
	static struct pscmp *DCFS_Af;
	int DCFS_Af_size;
	int NF_max=0, NP_max=0;
	static int time_in=0;
	static int **nn;	//nearest neighbours points.
	static double **interp_DCFS;
	static double *mycmb=NULL;
	static double **cmb_cumu;
	int n_withslimodel, n_withoutslimodel;
//	FILE *fout;
	time_in+=1;

	//----------------------------------------------------------------------//
	// 		 		This initial part will be executed only once			//
	//----------------------------------------------------------------------//


	if (time_in==1) {
		print_logfile("\nSetting up variables for calculating Coulomb stress fields...\n");
		print_screen("\nSetting up variables for calculating Coulomb stress fields...\n");

		//-----calculate neighbouring points--------//

		if (gridpoints_err==1) {
			nn=imatrix(1,NgridT,1,6);
			nearest_neighbours(NgridT, crst.nLat,crst.nLon,crst.nD, nn);
			interp_DCFS=dmatrix(1,NgridT,1,2);
			mycmb=dvector(1,NgridT);
		}


		switch (vary_recfault==1){
		case 1:
			//allocate memory to strike0, dip0, rake0.
			//if vary_recfault==1, receiver fault changes at each MC iteration.
			if (which_recfault==0) {
			//different focal mechanisms selected if MC iterations of focal mechanism should be used. However, if which_recfault!=0, a single mechanism should be used (focmec[X][which_recfault]).
				strike0=dvector(0,crst.nofmzones-1);
				dip0=dvector(0,crst.nofmzones-1);
				rake0=dvector(0,crst.nofmzones-1);
			}

			else {
				//in this case a single foc. mec. is needed.
				strike0=malloc(sizeof(double));
				dip0=malloc(sizeof(double));
				rake0=malloc(sizeof(double));
			}
			break;

		case 0:
			// if a fixed mechanism is used, the receiver fault from crst.str0, crst.dip0 should be used.
			// crst.str0[0] contains the regional mechanism, that should be used if no spatially variable mech. is given ((*crst).variable_fixmec=0).
			// crst.str0[1...NP] contain the foc. mech. for individual grid points, which should be used if (*crst).variable_fixmec=0.
			// This is achieved by passing crst.str0+fm_offset to the function "resolve_DCFS" later on.
			fm_offset= crst.variable_fixmec ? 1 : 0;
			break;

		default:
			break;
	}


		//-----------------------------------------------------------------//
		//					Afterslip (initial setup)					   //
		//-----------------------------------------------------------------//

		if (afterslip!=0){

			NTSeff=(multisnap)? NTScont : 1;
			DCFS_Af_size= NTSeff*NA;
			DCFS_Af= pscmp_array(0,DCFS_Af_size);

			a_ev=0;	//counter for DCFS_Af
			nfaults=0;	//counter for eqkfmAf

			if (multisnap) {
				cmb_cumu=dmatrix(0,NA-1,1,NgridT);
				for (int a=0; a<NA; a++){
					for (int n=1; n<=NgridT; n++) cmb_cumu[a][n]=0.0;
				}
			}

			//Find element of AllCoeff which should also be used for afterslip:
			AllCoeffaft=AllCoeff;
			for (int i=0; i<Nmain; i++){
				if (!AllCoeffaft->hasafterslip) {
					AllCoeffaft=AllCoeffaft->next;
				}

				else{

					for (int i=0; i<NTSeff; i++){
						DCFS_Af[a_ev+i].NF=AllCoeffaft->NF;
						DCFS_Af[a_ev+i].cmb= dvector(1,NgridT);	//only allocated stuff needed by OkadaCoeff2....
						DCFS_Af[a_ev+i].S=d3tensor(1,NgridT,1,3,1,3);
						DCFS_Af[a_ev+i].nsel=eqkfmAf[nfaults].nsel;
						DCFS_Af[a_ev+i].which_pts=eqkfmAf[nfaults].selpoints;
					}

					Coeffs_st=AllCoeffaft->Coeffs_st;
					Coeffs_dip=AllCoeffaft->Coeffs_dip;
					Coeffs_open=AllCoeffaft->Coeffs_open;

					for (int i=0; i<NTSeff; i++)	{
						for (int nf=0; nf<DCFS_Af[a_ev].NF; nf++){
							//assign slip values for each snapshot and fault into slip_X (this is needed because okadaCoeff2DCFS uses them):
							eqkfmAf[nfaults+nf].slip_str= (eqkfmAf[nfaults+nf].allslip_str)? eqkfmAf[nfaults+nf].allslip_str[i] : NULL;
							eqkfmAf[nfaults+nf].slip_dip= (eqkfmAf[nfaults+nf].allslip_dip)? eqkfmAf[nfaults+nf].allslip_dip[i] : NULL;
							eqkfmAf[nfaults+nf].open= (eqkfmAf[nfaults+nf].allslip_open)? eqkfmAf[nfaults+nf].allslip_open[i] : NULL;
						}
						//todo make this work for splines==1 too...
						okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, Coeffs_open, DCFS_Af[a_ev+i], eqkfmAf+nfaults); //todo free memory used by *AllCoeff;
						if (vary_recfault==0) {
							resolve_DCFS(DCFS_Af[a_ev+i], crst, crst.str0+fm_offset, crst.dip0+fm_offset, NULL, 1);
							//resolve_DCFS(DCFS_Af[a_ev+i], crst, crst.str0+fm_offset, crst.dip0+fm_offset, crst.rake0+fm_offset, 0);	//fixme one line or the other
							free_d3tensor(DCFS_Af[a_ev+i].S, 1,NgridT,1,3,1,3);
							if (multisnap && i<NTSeff){
								for (int n=1; n<=NgridT; n++) cmb_cumu[a_ev/NTSeff][n]+=DCFS_Af[a_ev+i].cmb[n];
							}
						}
					}
					nfaults+=DCFS_Af[a_ev].NF;	//counter for eqkfmAf
					a_ev+=NTSeff;	//counter for DCFS_Af
					AllCoeffaft=AllCoeffaft->next;

				}
			}
		}
		print_screen("done.\n");
	}

	//-----------------------------------------------------------------//
	//					Mainshock(initial setup)					   //
	//(this part executed every time a new set of slip models is used) //
	//-----------------------------------------------------------------//

	if (time_in==1 || refresh){
		if (afterslip!=2){
			NFsofar=0;
			temp=AllCoeff;
			for (int i=0; i<Nmain; i++){

				//Don't do anything if event is outside data time period:
				if (DCFS[temp->which_main].t <tdata0 || DCFS[temp->which_main].t>tdata1){
					NFsofar+=temp->NF;
					temp=temp->next;
					continue;
				}

				if (eqkfm0[NFsofar].is_slipmodel){
					//calculate stress tensor at each grid point:
					Coeffs_st=temp->Coeffs_st;
					Coeffs_dip=temp->Coeffs_dip;
					Coeffs_open=temp->Coeffs_open;
					okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, Coeffs_open, DCFS[temp->which_main], eqkfm0+NFsofar);

					//resolve coefficients if receiver faults don't change between iterations:
					switch (vary_recfault){
						case 0:
							resolve_DCFS(DCFS[temp->which_main], crst, crst.str0+fm_offset, crst.dip0+fm_offset, NULL, 1);	//fixme choose one
							//resolve_DCFS(DCFS[temp->which_main], crst, crst.str0+fm_offset, crst.dip0+fm_offset, crst.rake0+fm_offset, 0);
							if (gridpoints_err){
								int eq1=temp->which_main;

								//need to copy field into cmb0, so that each iteration will smooth the original field (and not the one smoothed in prev. iteration).
								for (int i=1; i<=NgridT; i++) mycmb[i]=0.0;
								for (int i=1; i<=DCFS[eq1].nsel; i++) mycmb[DCFS[eq1].which_pts[i]]=DCFS[eq1].cmb[i];
								interp_nn(NgridT,crst.nLat, crst.nLon, crst.nD, mycmb,interp_DCFS,0,nn);
								DCFS[eq1].cmb0=dvector(1,DCFS[eq1].nsel);
								DCFS[eq1].Dcmb=dvector(1,DCFS[eq1].nsel);
								for (int i=1; i<=DCFS[eq1].nsel; i++){
									DCFS[eq1].cmb0[i]=0.5*(interp_DCFS[DCFS[eq1].which_pts[i]][1]+interp_DCFS[DCFS[eq1].which_pts[i]][2]);
									DCFS[eq1].Dcmb[i]=fabs(interp_DCFS[DCFS[eq1].which_pts[i]][1]-interp_DCFS[DCFS[eq1].which_pts[i]][2]);
								}
							}

							break;
						case 2:
							// todo [coverage] this block is never tested
							DCFScmbopt(DCFS, temp->which_main, crst);	//NB this does not take into account stress from afterslip, assuming that from mainshocks is much larger this is ok.
							break;
						default:
							break;
					}

				}

				else{
					//prepare isotropic stress fields:
					if (!eqkfm0[NFsofar].is_slipmodel){
						int eq1=temp->which_main;
						isoDCFS(DCFS[eq1], eqkfm0[NFsofar]);
						if (gridpoints_err){
							for (int i=1; i<=NgridT; i++) mycmb[i]=0.0;
							for (int i=1; i<=DCFS[eq1].nsel; i++) mycmb[DCFS[eq1].which_pts[i]]=DCFS[eq1].cmb[i];
							interp_nn(NgridT,crst.nLat, crst.nLon, crst.nD, mycmb,interp_DCFS,0,nn);
							DCFS[eq1].cmb0=dvector(1,DCFS[eq1].nsel);
							DCFS[eq1].Dcmb=dvector(1,DCFS[eq1].nsel);
							for (int i=1; i<=DCFS[eq1].nsel; i++){
								DCFS[eq1].cmb0[i]=0.5*(interp_DCFS[DCFS[eq1].which_pts[i]][1]+interp_DCFS[DCFS[eq1].which_pts[i]][2]);
								DCFS[eq1].Dcmb[i]=fabs(interp_DCFS[DCFS[eq1].which_pts[i]][1]-interp_DCFS[DCFS[eq1].which_pts[i]][2]);
							}
						}
					}
				}
				NFsofar+=temp->NF;
				temp=temp->next;
			}
		}
	}


	//----------------------------------------------------------------------//
	//	 		From here on will be executed at every iteration			//
	//----------------------------------------------------------------------//


	//At the moment afterslip and OOPs at the same time are not implemented.
	if (afterslip==1 && vary_recfault==2) {
		print_screen("*Error: function calculateDCFSperturbed doesn't know how to calculate OOPS when afterslip is included!!*\n");
		print_logfile("*Error: function calculateDCFSperturbed doesn't know how to calculate OOPS when afterslip is included!!*\n");
		return;
	}

	//pick receiver faults from catalog of focal mechanisms:
	if (vary_recfault==1){
		//if vary_recfault==1, receiver fault changes at each MC iteration.
		if (which_recfault==0) {
			//randomly pick a mechanism for each zone:
			for (int fmzone=0; fmzone<crst.nofmzones; fmzone++){
				first=fmzoneslim[fmzone];
				last=fmzoneslim[fmzone+1]-1;
				rand= (int) ((last-first)*ran1(seed)+first);
				*seed=-*seed;
				strike0[fmzone]=focmec[1][rand];
				dip0[fmzone]=focmec[2][rand];
				rake0[fmzone]=focmec[3][rand];	//only used for splines==1 (see below).
			}
		}
		else {
			// use the foc. mec. for this iteration (this is done when all focal mechanisms should be sampled; only activated in main.c if nofmzones=1).
			// todo [coverage] this block is never tested
			*strike0=focmec[1][which_recfault];
			*dip0=focmec[2][which_recfault];
			*rake0=focmec[3][which_recfault];	//only used for splines==1 (see below).
		}
	}

	//------------------------------------------//
	//	calculated stress field from afterslip: //
	//------------------------------------------//

	if (afterslip==0){
		for (int l=0; l<NTScont; l++){
			if (times[l] <tdata0 || times[l]>tdata1) continue;
		// todo [coverage] this line is never tested
			for (int n=1; n<=NgridT; n++) if (DCFSrand) DCFSrand[l][n]=0.0;	//todo should DCFSrand just set to NULL if afterslip==0?
		}
	}
	else {
		i=0;
		for (int a=0; a<NA; a++){

			if (multisnap==0){

				if (vary_recfault==1) resolve_DCFS(DCFS_Af[a], crst, strike0, dip0, NULL, 1);
				//if (vary_recfault==1) resolve_DCFS(DCFS_Af[a], crst, strike0, dip0, rake0, 0);	//fixme choose a line
				for (int l=0; l<NTScont; l++) {
					if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;

					//loop over events with afterslip:
					for (int n=1; n<=NgridT; n++) DCFSrand[l][n]=eqkfmAf[i].tevol[l]*DCFS_Af[a].cmb[n];	//FIXME not eqkfmAf[0].
				}
				i+=DCFS_Af[a].NF;
			}

			else{
				if (vary_recfault==1){
					//fixme cmb_cumu[0];
					for (int n=1; n<=NgridT; n++) cmb_cumu[a][n]=0.0;
					for (int l=0; l<NTScont; l++) {
						resolve_DCFS(DCFS_Af[a*NTScont+l], crst, strike0, dip0, NULL, 1); //fixme choose one
						//resolve_DCFS(DCFS_Af[a*NTScont+l], crst, strike0, dip0, rake0, 0);
						if (l<NTScont-1) for (int n=1; n<=NgridT; n++) cmb_cumu[a][n]+=DCFS_Af[a*NTScont+l].cmb[n];
					}
				}

				for (int l=0; l<NTScont; l++) {
					if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
					for (int n=1; n<=NgridT; n++) {
						DCFSrand[l][n]= (fabs(cmb_cumu[a][n])>DCFS_cap) ? (DCFS_cap/fabs(cmb_cumu[a][n]))*DCFS_Af[a*NTScont+l].cmb[n] : DCFS_Af[a*NTScont+l].cmb[n];
					}
				}
			}
		}
	}

	//------------------------------------------------------------------------------------------------------//
	//	calculated stress field from mainshocks (i.e. events for which a non trivial slip model is used):   //
	//------------------------------------------------------------------------------------------------------//

	temp=AllCoeff;
	NFsofar=0;
	for (int i=0; i<Nmain; i++){

		if (DCFS[temp->which_main].t <tdata0 || DCFS[temp->which_main].t>tdata1){
			NFsofar+=temp->NF;
			temp=temp->next;
			continue;
		}

		else{
			if (afterslip==2){
				for (int n=1; n<=NgridT; n++) DCFS[temp->which_main].cmb[n]=0.0;	//don't need to fill in tensor.
			temp=temp->next;
			}
			else{
				if (eqkfm0[NFsofar].is_slipmodel){
					if (vary_recfault==1){	// In this case stress tensor has to be resolved on a different plane at each iteration:
						resolve_DCFS(DCFS[temp->which_main], crst, strike0, dip0, NULL, 1); //fixme choose a line
						//resolve_DCFS(DCFS[temp->which_main], crst, strike0, dip0, rake0, 0);
					}
					if (gridpoints_err==1) {
						//if vary_recfault==0, the field hasn't changed and cmb0 should be used.
						smoothen_DCFS(DCFS[temp->which_main], crst.nLat, crst.nLon, crst.nD, seed, vary_recfault==0, nn);
					}
				}
				else {
					//if the earthquake does not have a slip model, the isotropic field does not change between iterations and has already been calculated.
					//only need to calculated error from gridpoints.
					if (gridpoints_err==1) smoothen_DCFS(DCFS[temp->which_main], crst.nLat, crst.nLon, crst.nD, seed, 1, nn);
				}

				NFsofar+=temp->NF;
				temp=temp->next;
			}
		}
	}
}

//todo move soewhere else
void smoothen_DCFS(struct pscmp DCFS, int nlat, int nlon, int nd, long *seed, int use_cmb0, int **nn){
	//can use DCFS.cmb (which gets overwritten) or DCFS.cmb0 (which is preserved) -> latter is useful if grid point smoothing is only source of uncertainty.
	//with this option, DCFS.cmb0= mean value of field; DCFS.cmb=range of values.

	double randcmb;
	double *mycmb, **interp_DCFS;
	int NgridT=nlat*nlon*nd;

	if (!use_cmb0){
		mycmb=dvector(1,NgridT);
		interp_DCFS=dmatrix(1,NgridT,1,2);
		for (int i=1; i<=NgridT; i++) mycmb[i]=0.0;
		for (int i=1; i<=DCFS.nsel; i++) mycmb[DCFS.which_pts[i]]=DCFS.cmb[i];

		interp_nn(NgridT, nlat, nlon, nd, mycmb, interp_DCFS, 0, nn);

		for (int i=1; i<=DCFS.nsel; i++) {
			int p=DCFS.which_pts[i];
			randcmb=interp_DCFS[p][1]+ran1(seed)*(interp_DCFS[p][2]-interp_DCFS[p][1]);
			*seed=-*seed;
			DCFS.cmb[i]=randcmb;
		}

		free_dvector(mycmb,1,NgridT);
		free_dmatrix(interp_DCFS,1,NgridT,1,2);
	}
	else{
		for (int i=1; i<=DCFS.nsel; i++) {
			randcmb=-0.5*DCFS.Dcmb[i]+ran1(seed)*DCFS.Dcmb[i];
			*seed=-*seed;
			DCFS.cmb[i]=DCFS.cmb0[i]+randcmb;
		}
	}
}

//void smoothen_vector(int NgridT, int nLat, int nLon, int nD, double *values, long *seed, int **nn, int return_range){
//	double randcmb;
//	double **interp_DCFS=dmatrix(1,NgridT, 1, 2);
//
//	interp_nn(NgridT,nLat,nLon,nD,values,interp_DCFS,0,nn);
//	for (int i=1; i<=NgridT; i++) {
//		if (return_range) values[i]=interp_DCFS[i][2]-interp_DCFS[i][1];
//		else {
//			randcmb=interp_DCFS[i][1]+ran1(seed)*(interp_DCFS[i][2]-interp_DCFS[i][1]);
//			*seed=-*seed;
//			values[i]=randcmb;
//		}
//	}
//	free_dmatrix(interp_DCFS,1,NgridT, 1, 2);
//}
