/*
 * calculateDCFSperturbed.c
 *
 *  Created on: Mar 19, 2013
 *      Author: camcat
 */

#include "calculateDCFSperturbed.h"

void calculateDCFSperturbed(double **DCFSrand, struct pscmp *DCFS, struct eqkfm *eqkfmAf,
							struct eqkfm *eqkfm0, struct eqkfm *eqkfm1, struct flags flag,
							double *tevol, double *times, int Nmain, struct crust crst,
							struct Coeff_LinkList *AllCoeff, int NTScont, int NTSdisc,
							double **focmec, int *fmzoneslim, int NFM, long *seed,
							int *fm_number, double tdata0, double tdata1, double H,
							int refresh, int which_recfault) {

/* DCFSrand[i][j] contains the ith stress change at gridpoint j due to continuous processes (modeled linearly between time steps).
 * DCFS[k].cmp[j] contains the stress change due to kth event at gridpoint j (modeled as a step).
 * NTScont is the total number of time steps for continuos process;
 * NTSdisc is the total number of time steps for step like process (i.e. number of triggering earthquakes).
 * focmec contains sample of focal mechanisms, NFM is its length (1->NFM)
 * afterslip and aftershocks are flags indicating if these processes should be included.
 * vary_recfault, vary_slipmodel are flags indicating which sources of uncertainties should be included. 0: use fix planes; 1: vary foc. mec. (use which_fm, or if which_fm==0, choose random one); 2: use OOPs.
 * which_recfault gives index of foc.mec to select; if set to 0, choose random one (useful if NFM>>nsur). if vary_recfault=0, this is meaningless.
 *
 * NB: not general, as it assumes mainshock is the first event. todo check/fix this.
 */

	//flags:
	int	afterslip=flag.afterslip, \
		aftershocks=flag.aftershocks, \
		vary_recfault=flag.err_recfault, \
		vary_slipmodel=flag.err_slipmodel, \
		new_slipmodel=flag.new_slipmodel, \
		gridpoints_err=flag.err_gridpoints, \
		splines=flag.splines, \
		afterslip_errors=flag.err_afterslipmodel,	\
		full_field=(flag.full_field==2), \
		aftershocks_fixedmec=flag.aftershocks_fixedmec;

	static double *strike0, *dip0, *rake0;
	double slip;
	double ***Stemp;
	static double lat0, lon0;
	int rand;
	int NgridT=crst.N_allP;
	int NFsofar=0;
	int last, first;
	static struct eqkfm *eqkfm2;
	static struct eqkfm *eqkfm_noise;
	static struct eqkfm *eqkfm2A;
	struct Coeff_LinkList *temp;
	float ***Coeffs_st, ***Coeffs_dip;	//coefficients for tensor;
	static float **Coeff_ResS, **Coeff_ResD;	//coefficients for cmb (resolved).
	static struct pscmp *DCFS_Af, DCFS_Af_noise;
	int DCFS_Af_size;
	int NF_max=0, NP_max=0;
	static int time_in=0;
	static int **nn;	//nearest neighbours points.
	static double **interp_DCFS;
	static double *mycmb=NULL;
	static double *cmb_cumu;
	int n_withslimodel, n_withoutslimodel;
//	FILE *fout;
	time_in+=1;

	if (time_in==1) {

		if (flog) fprintf(flog, "\nSetting up variables for calculating perturbed Coulomb fields.\n");

		if (vary_slipmodel && afterslip_errors && splines && gridpoints_err){
			if (verbose_level>0) printf("** Warning: flags vary_slipmodel, afterslip_errors, splines, gridpoints_err set to 1 -> slow! (calculateDCFSperturbed.c)**\n");
		}

		new_slipmodel=1;
		//-----calculate neighbouring points--------//

		if (gridpoints_err==1) {
			nn=imatrix(1,NgridT,1,6);
			nearest_neighbours(NgridT, crst.nLat,crst.nLon,crst.nD, nn);
			interp_DCFS=dmatrix(1,NgridT,1,2);
			mycmb=dvector(1,NgridT);
		}

		if (vary_recfault==1 && which_recfault==0){
			strike0=dvector(0,crst.nofmzones-1);
			dip0=dvector(0,crst.nofmzones-1);
			rake0=dvector(0,crst.nofmzones-1);
		}
		else {
			strike0=malloc(sizeof(double));
			dip0=malloc(sizeof(double));
			rake0=malloc(sizeof(double));
		}

		//-----------------------------------------------------------------//
		//							Afterslip							   //
		//-----------------------------------------------------------------//

		if (afterslip!=0){
		//fixme afterslip should have its own structure.

			if (splines) {
				cmb_cumu=dvector(1,NgridT);
				for (int n=1; n<=NgridT; n++) cmb_cumu[n]=0.0;
			}

			DCFS_Af_size= (splines)? NTScont+1: 1;
			DCFS_Af= pscmp_array(0,DCFS_Af_size-1);

			DCFS_Af[0].NF=AllCoeff->NF;
			DCFS_Af[0].cmb=dvector(1,NgridT);	//only allocated stuff needed by OkadaCoeff2....
			DCFS_Af[0].S=d3tensor(1,NgridT,1,3,1,3);
			DCFS_Af[0].nsel=eqkfmAf[0].nsel;
			DCFS_Af[0].which_pts=eqkfmAf[0].selpoints;

			for (int i=1; i<DCFS_Af_size; i++){
				DCFS_Af[i].NF=DCFS_Af[0].NF;
				DCFS_Af[i].cmb= dvector(1,NgridT);
				DCFS_Af[i].S=d3tensor(1,NgridT,1,3,1,3);
				DCFS_Af[i].nsel=eqkfmAf[0].nsel;
				DCFS_Af[i].which_pts=eqkfmAf[0].selpoints;
			}

			if (splines && afterslip_errors && vary_slipmodel && vary_recfault==0){
				okadaCoeff_resolve(AllCoeff[0], &Coeff_ResS, &Coeff_ResD, crst, &crst.str0, &crst.dip0, &crst.rake0);	//to reduce calculations.
			}

			if (vary_slipmodel && afterslip_errors){
				for (int f=0; f<AllCoeff->NF; f++) NP_max=fmax(NP_max,eqkfmAf[f].np_di*eqkfmAf[f].np_st);
				if (splines){
					eqkfm_noise=eqkfm_array(0, AllCoeff->NF);
					for (int f=0; f<AllCoeff->NF; f++) {
						eqkfm_noise[f].slip_dip=dvector(1,NP_max);
						eqkfm_noise[f].slip_str=dvector(1,NP_max);
					}
					DCFS_Af_noise.NF=AllCoeff->NF;
					DCFS_Af_noise.cmb=dvector(1,NgridT);	//only allocated stuff needed by OkadaCoeff2....
					DCFS_Af_noise.S=d3tensor(1,NgridT,1,3,1,3);
					DCFS_Af_noise.nsel=eqkfmAf[0].nsel;
					DCFS_Af_noise.which_pts=eqkfmAf[0].selpoints;
				}
				else {
					eqkfm2A=eqkfm_array(0, AllCoeff->NF);		//todo do not assume first event is the one from which afterslip geometry is taken.
					for (int f=0; f<AllCoeff->NF; f++) {
						eqkfm2A[f].slip_dip=dvector(1,NP_max);
						eqkfm2A[f].slip_str=dvector(1,NP_max);
					}
				}
			}

			else{
				Coeffs_st=AllCoeff->Coeffs_st;
				Coeffs_dip=AllCoeff->Coeffs_dip;
				for (int i=0; i<DCFS_Af_size; i++)	{
					okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, DCFS_Af[i], eqkfmAf+i*DCFS_Af[0].NF, crst, NULL, NULL, 1); //todo free memory used by *AllCoeff; todo make this work for splines==1 too...
					if (vary_recfault==0) {
						resolve_DCFS(DCFS_Af[i], crst, &crst.str0, &crst.dip0, NULL, 1);
						if (splines && i<DCFS_Af_size-1){
							for (int n=1; n<=NgridT; n++) cmb_cumu[n]+=DCFS_Af[i].cmb[n];
						}
					}
				}
				if (gridpoints_err && afterslip_errors){
					DCFS_Af[0].cmb0=dvector(1,NgridT);	//only allocated stuff needed by OkadaCoeff2....
					DCFS_Af[0].Dcmb=dvector(1,NgridT);	//only allocated stuff needed by OkadaCoeff2....
					if (!splines && vary_recfault==0){
						for (int i=1; i<=NgridT; i++) DCFS_Af[0].cmb0[i]=DCFS_Af[0].cmb[i];
						smoothen_vector(NgridT, crst.nLat, crst.nLon, crst.nD, DCFS_Af[0].cmb, seed, nn, 1);
						for (int i=1; i<=NgridT; i++) DCFS_Af[0].Dcmb[i]=DCFS_Af[0].cmb[i];	// contains range.
					}
					if (splines) {
						for (int i=0; i<DCFS_Af_size; i++){
							DCFS_Af[i].cmb0=DCFS_Af[0].cmb0;
							DCFS_Af[i].Dcmb=DCFS_Af[0].Dcmb;
						}
					}
				}
			}
		}

		//-----------------------------------------------------------------//
		//							Aftershocks							   //
		//-----------------------------------------------------------------//

		if (aftershocks){
			lat0=crst.lat0;
			lon0=crst.lon0;
			n_withslimodel=n_withoutslimodel=0;
			for (int eq1=0; eq1<NTSdisc; eq1++){
				if (!eqkfm1[eq1].is_mainshock && eqkfm1[eq1].nsel!=0){
					if (eqkfm1[eq1].is_slipmodel){
						n_withslimodel+=1;
						okadaDCFS(DCFS[eq1], eqkfm1+eq1, 1, crst, NULL, NULL, 1);
						if (eqkfm1[eq1].whichfm==0){	//also calculate S tensor for second slip model.
							DCFS[eq1].S1=DCFS[eq1].S;
							DCFS[eq1].S=d3tensor(1, DCFS[eq1].nsel, 1,3,1,3);
							eqkfm1[eq1].whichfm=2;
							double temp_strslip=eqkfm1[eq1].slip_str[1];
							double temp_dipslip=eqkfm1[eq1].slip_dip[1];
							eqkfm1[eq1].slip_str[1]=eqkfm1[eq1].slip_str[2];
							eqkfm1[eq1].slip_dip[1]=eqkfm1[eq1].slip_dip[2];
							okadaDCFS(DCFS[eq1], eqkfm1+eq1, 1, crst, NULL, NULL, 1);
							eqkfm1[eq1].slip_str[1]=temp_strslip;
							eqkfm1[eq1].slip_dip[1]=temp_dipslip;
						}
					}

					else {
						n_withoutslimodel+=1;
						if (full_field){
							if (aftershocks_fixedmec){
								eqkfm1[eq1].str1=crst.str0;
								eqkfm1[eq1].dip1=crst.dip0;
								eqkfm1[eq1].rake1=crst.rake0;
								eqkfm1[eq1].whichfm=1;
								WellsCoppersmith(eqkfm1[eq1].mag, eqkfm1[eq1].rake1, &(eqkfm1[eq1].L), &(eqkfm1[eq1].W), &slip);
								slip=pow(10.0,1.5*(eqkfm1[eq1].mag+6.0))/(pow(10,12)*crst.mu*eqkfm1[eq1].W*eqkfm1[eq1].L);
								eqkfm1[eq1].slip_str=dvector(1,1);
								eqkfm1[eq1].slip_dip=dvector(1,1);
								eqkfm1[eq1].slip_str[1]=slip*cos(DEG2RAD*eqkfm1[eq1].rake1);
								eqkfm1[eq1].slip_dip[1]=slip*sin(DEG2RAD*eqkfm1[eq1].rake1);
								okadaDCFS(DCFS[eq1], eqkfm1+eq1, 1, crst, NULL, NULL, 1);
							}
							else {
								eqkfm1[eq1].slip_str=dvector(1,1);
								eqkfm1[eq1].slip_dip=dvector(1,1);
								eqkfm1[eq1].pos_s=dvector(1,1);
								eqkfm1[eq1].pos_d=dvector(1,1);
							}
						}

					else {
							isoDCFS(DCFS[eq1], eqkfm1[eq1]);
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
				}
			}

			if (flog){
				fprintf(flog,"%d events have known focal mechanism, which will be used.\n", n_withslimodel);
				fprintf(flog,"%d events do not have a known focal mechanism. ", n_withoutslimodel);
				if (n_withoutslimodel){
					if (full_field) {
						fprintf(flog,"will use ");
						if (aftershocks_fixedmec) fprintf(flog,"fixed mechanism (strike=%.2lf, dip=%.2lf).\n", crst.str0, crst.dip0);
						else fprintf(flog,"Monte Carlo sampling from catalog of focal mechanisms.\n");
					}
					else fprintf(flog,"will use isotropic field.\n");
				}
				fflush(flog);
			}
		}

	}
		//-----------------------------------------------------------------//
		//							Mainshock							   //
		//-----------------------------------------------------------------//

	if (time_in==1 || refresh){
		if (afterslip!=2){
			if (vary_slipmodel){
				//find largest no. of faults and largest no of patches needed by eqkfm2.
				temp=AllCoeff;
				for (int i=0; i<Nmain; i++) {
					NF_max=fmax(NF_max,temp->NF);
					for (int f=0; f<temp->NF; f++) NP_max=fmax(NP_max,eqkfm0[f+NFsofar].np_di*eqkfm0[f+NFsofar].np_st);
					NFsofar+=temp->NF;
					temp=temp->next;
				}
				if (time_in!=1) {
					for (int f=0; f<NF_max; f++) {
						free_dvector(eqkfm2[f].slip_dip,1,0);
						free_dvector(eqkfm2[f].slip_str,1,0);
					}
					free_eqkfmarray(eqkfm2, 0, 0);
				}
				eqkfm2=eqkfm_array(0, NF_max);		//todo also initialize vectors inside (with size of the largest), if they are to be recycled...
				for (int f=0; f<NF_max; f++) {
					eqkfm2[f].slip_dip=dvector(1,NP_max);
					eqkfm2[f].slip_str=dvector(1,NP_max);
				}
			}
			else{
				NFsofar=0;
				temp=AllCoeff;
				for (int i=0; i<Nmain; i++){
					if (eqkfm0[NFsofar].is_slipmodel){
						Coeffs_st=temp->Coeffs_st;
						Coeffs_dip=temp->Coeffs_dip;
						okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, DCFS[temp->which_main], eqkfm0+NFsofar, crst, NULL, NULL, 1);
					}
					NFsofar+=temp->NF;
					temp=temp->next;
				}
			}

			//prepare isotropic fields:
			NFsofar=0;
			temp=AllCoeff;
			for (int i=0; i<Nmain; i++){
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
				NFsofar+=temp->NF;
				temp=temp->next;
			}
		}
	}


	if (afterslip==1 && vary_recfault==2) {
		printf("*Error: function calculateDCFSperturbed doesn't know how to calculate OOPS when afterslip is included!!*\n");
		return;
	}

	switch (vary_recfault){
		case 1:
			if (which_recfault==0) {
				for (int fmzone=0; fmzone<crst.nofmzones; fmzone++){
					first=fmzoneslim[fmzone];
					last=fmzoneslim[fmzone+1]-1;
					rand= (int) ((last-first)*ran1(seed)+first); ///*********************************
					*seed=-*seed;
					strike0[fmzone]=focmec[1][rand];
					dip0[fmzone]=focmec[2][rand];
					rake0[fmzone]=focmec[3][rand];	//only used for splines==1 (see below).
				}
				if (fm_number){
					if (crst.nofmzones==1) *fm_number=rand;
					else {
						*fm_number=0;
						if (verbose_level>1) printf("**Warning: multiple focal mechanisms catalogs available, can not print out single value! (calculateDCFSperturbed.c).**\n");
						if (flog) fprintf(flog,"**Warning: multiple focal mechanisms catalogs available, can not print out single value! (calculateDCFSperturbed.c).**\n");
					}
				}
			}
			else {
				*strike0=focmec[1][which_recfault];
				*dip0=focmec[2][which_recfault];
				*rake0=focmec[3][which_recfault];	//only used for splines==1 (see below).
				if (fm_number) *fm_number=which_recfault;
			}
			break;

		case 0:
			*strike0= crst.str0;
			*dip0=crst.dip0;
			*rake0=crst.rake0;	//only used for splines==1 (see below).
			if (fm_number!=0) *fm_number=0;
			break;

		case 2:
			break;	//receiver fault varies for each point (OOP).
	}

	//------------------------------------------//
	//	calculated stress field from afterslip: //
	//------------------------------------------//

	if (afterslip==0){
		for (int l=0; l<NTScont; l++){
			if (times[l] <tdata0 || times[l]>tdata1) continue;
			for (int n=1; n<=NgridT; n++) if (DCFSrand) DCFSrand[l][n]=0.0;
		}
	}
	else {
		if (times[0]<tdata1 && times[NTScont-1]>=tdata0){

			if (splines==0){
				if (vary_slipmodel && afterslip_errors){
					if (new_slipmodel) {
						for (int f=0; f<DCFS_Af[0].NF; f++) {
							suomod1_hf(eqkfmAf[f], eqkfm2A+f, H, seed, 0);
						}
					}
					//todo create a AllCoeff structure for afterslip?
					//todo rake should be optimized considering *all* prev. stress...
					//todo deal with case vary_recfault==2?
					okadaCoeff2DCFS(AllCoeff->Coeffs_st, AllCoeff->Coeffs_dip, DCFS_Af[0], eqkfm2A, crst, strike0, dip0, 0);
				}
				else {
					if (vary_recfault!=0) resolve_DCFS(DCFS_Af[0], crst, strike0, dip0, NULL, 1);
				}

				if (gridpoints_err && afterslip_errors) smoothen_DCFS(DCFS_Af[0], crst.nLat, crst.nLon, crst.nD, seed, (!vary_recfault && !vary_slipmodel), nn);
				for (int l=0; l<NTScont; l++) {
					if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
					for (int n=1; n<=NgridT; n++) DCFSrand[l][n]=tevol[l]*DCFS_Af[0].cmb[n];
				}
			}

			else{
				if (afterslip_errors && vary_slipmodel){
					if (vary_recfault!=0) okadaCoeff_resolve(AllCoeff[0], &Coeff_ResS, &Coeff_ResD, crst, strike0, dip0, rake0);	//to reduce calculations.
					if (new_slipmodel){
						for (int f=0; f<DCFS_Af[0].NF; f++) {
							suomod1_hf(eqkfmAf[DCFS_Af[0].NF*NTScont+f], eqkfm_noise+f, H, seed, 1);	//calculate noise from cumulative slip.
							//copy_eqkfm_all(eqkfmAf[DCFS_Af[0].NF*NTScont+f], eqkfm_noise+f);
						}
						okadaCoeff2DCFS(AllCoeff->Coeffs_st, AllCoeff->Coeffs_dip, DCFS_Af_noise, eqkfm_noise, crst, strike0, dip0, 0);
						if (gridpoints_err){
							resolvedCoeff2DCFS(Coeff_ResS, Coeff_ResD, DCFS_Af[0], eqkfmAf+DCFS_Af[0].NF*NTScont, crst);
							for (int n=1; n<=NgridT; n++) cmb_cumu[n]=DCFS_Af[0].cmb[n]+DCFS_Af_noise.cmb[n];
							smoothen_vector(NgridT,crst.nLat, crst.nLon, crst.nD,cmb_cumu,seed, nn, 1);	//gives range of possible values.
						}
					}
					for (int n=1; n<=NgridT; n++) cmb_cumu[n]=0.0;
					for (int l=0; l<NTScont; l++) {
						if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
						resolvedCoeff2DCFS(Coeff_ResS, Coeff_ResD, DCFS_Af[l], eqkfmAf+DCFS_Af[l].NF*l, crst);
						if (l<NTScont-1) {
							for (int n=1; n<=NgridT; n++) {
								cmb_cumu[n]+=DCFS_Af[l].cmb[n]+(eqkfmAf[DCFS_Af[l].NF*l].tot_slip/eqkfm_noise[0].tot_slip)*DCFS_Af_noise.cmb[n];
							}
						}
					}
					for (int l=0; l<NTScont; l++) {
						if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
						if (gridpoints_err){
							for (int n=1; n<=NgridT; n++){
								DCFS_Af[l].cmb0[n]=DCFS_Af[l].cmb[n]+(eqkfmAf[DCFS_Af[l].NF*l].tot_slip/eqkfm_noise[0].tot_slip)*DCFS_Af_noise.cmb[n];
								DCFS_Af[l].Dcmb[n]=cmb_cumu[n]*(eqkfmAf[DCFS_Af[l].NF*l].tot_slip/eqkfm_noise[0].tot_slip);		//range scaled by total slip.
							}
							smoothen_DCFS(DCFS_Af[l], crst.nLat, crst.nLon, crst.nD, seed, 1, nn);
							for (int n=1; n<=NgridT; n++) {
								DCFSrand[l][n]= (fabs(cmb_cumu[n])>DCFS_cap) ? (DCFS_cap/fabs(cmb_cumu[n]))*DCFS_Af[l].cmb[n] : DCFS_Af[l].cmb[n];
							}
						}
						else {
							for (int n=1; n<=NgridT; n++) {
								DCFSrand[l][n]=DCFS_Af[l].cmb[n]+(eqkfmAf[DCFS_Af[l].NF*l].tot_slip/eqkfm_noise[0].tot_slip)*DCFS_Af_noise.cmb[n];
								if (fabs(cmb_cumu[n])>DCFS_cap) {
									DCFSrand[l][n]=DCFSrand[l][n]*(DCFS_cap/fabs(cmb_cumu[n]));
								}
							}
						}
					}
				}

				else {
					if (vary_recfault!=0){
						for (int n=1; n<=NgridT; n++) cmb_cumu[n]=0.0;
						for (int l=0; l<NTScont; l++) {
							resolve_DCFS(DCFS_Af[l], crst, strike0, dip0, NULL, 1);
							if (l<NTScont-1) for (int n=1; n<=NgridT; n++) cmb_cumu[n]+=DCFS_Af[l].cmb[n];
						}
					}

					for (int l=0; l<NTScont; l++) {
						if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
						if (afterslip_errors && gridpoints_err){
							for (int n=1; n<=NgridT; n++){
								DCFS_Af[l].cmb0[n]=DCFS_Af[l].cmb[n];
								DCFS_Af[l].Dcmb[n]=cmb_cumu[n]*(eqkfmAf[DCFS_Af[0].NF*l].tot_slip/eqkfmAf[DCFS_Af[0].NF*NTScont].tot_slip);		//range scaled by total slip.
							}
							smoothen_DCFS(DCFS_Af[l], crst.nLat, crst.nLon, crst.nD, seed, 1, nn);
						}
						for (int n=1; n<=NgridT; n++) {
							DCFSrand[l][n]= (fabs(cmb_cumu[n])>DCFS_cap) ? (DCFS_cap/fabs(cmb_cumu[n]))*DCFS_Af[l].cmb[n] : DCFS_Af[l].cmb[n];
						}
					}
				}
			}
		}
	}

	//------------------------------------------------------------------------------//
	//					calculated stress field from aftershocks					//
	//------------------------------------------------------------------------------//

	if (aftershocks==1){

		for (int eq1=0; eq1<NTSdisc; eq1++){

			if (eqkfm1[eq1].t <tdata0 || eqkfm1[eq1].t>tdata1) continue;

			if (!eqkfm1[eq1].is_mainshock && eqkfm1[eq1].nsel!=0){
				if (eqkfm1[eq1].is_slipmodel) {
					if (eqkfm1[eq1].whichfm==0){
						//half of the time, swap them.
						if ((int) round(ran1(seed))){
							Stemp=DCFS[eq1].S;
							DCFS[eq1].S=DCFS[eq1].S1;
							DCFS[eq1].S1=Stemp;
						}
						*seed=-*seed;
					}
					resolve_DCFS(DCFS[eq1], crst, strike0, dip0, NULL, 1);
					if (gridpoints_err) smoothen_DCFS(DCFS[eq1], crst.nLat, crst.nLon, crst.nD, seed,0, nn);
				}
				else {
					if (full_field){
						if (!aftershocks_fixedmec){
							rand=(int) ((NFM-1)*ran1(seed)+1);	//fixme: should choose from correct area if multiple foc mec are available (or pick the closest?)
							*seed=-*seed;
							eqkfm1[eq1].str1=focmec[1][rand];
							eqkfm1[eq1].dip1=focmec[2][rand];
							eqkfm1[eq1].rake1=fmod(focmec[3][rand]+360.0,360.0);

							WellsCoppersmith(eqkfm1[eq1].mag, eqkfm1[eq1].rake1, &(eqkfm1[eq1].L), &(eqkfm1[eq1].W), &slip);
							slip=pow(10.0,1.5*(eqkfm1[eq1].mag+6.0))/(pow(10,12)*crst.mu*eqkfm1[eq1].W*eqkfm1[eq1].L);
							eqkfm1[eq1].np_di=eqkfm1[eq1].np_st=1;
							eqkfm1[eq1].pos_d[1]=eqkfm1[eq1].pos_s[1]=0.0;
							eqkfm1[eq1].slip_str[1]=slip*cos(DEG2RAD*eqkfm1[eq1].rake1);
							eqkfm1[eq1].slip_dip[1]=slip*sin(DEG2RAD*eqkfm1[eq1].rake1);
							okadaDCFS(DCFS[eq1], eqkfm1+eq1, 1, crst, strike0, dip0, 1);
						}

						if (vary_recfault!=2) resolve_DCFS(DCFS[eq1], crst, strike0, dip0, NULL, 1);
						else DCFScmbopt(DCFS, eq1, crst);

						if (gridpoints_err==1) smoothen_DCFS(DCFS[eq1], crst.nLat, crst.nLon, crst.nD, seed,0, nn);
					}
					else {
						if (gridpoints_err) smoothen_DCFS(DCFS[eq1], crst.nLat, crst.nLon, crst.nD, seed,1, nn);
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
			}
			else{
				if (eqkfm0[NFsofar].is_slipmodel){
					Coeffs_st=temp->Coeffs_st;
					Coeffs_dip=temp->Coeffs_dip;
					if (vary_slipmodel){
						if (new_slipmodel){
							for (int f=0; f<temp->NF; f++) {
								//add high freq. noise.
								//if (time_in!=1) suomod1_cleanup(eqkfm2+f);
								suomod1_hf(eqkfm0[f+NFsofar], eqkfm2+f, H, seed, 0);
							}
						}
						if (vary_recfault==2) {
							okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, DCFS[temp->which_main], eqkfm2, crst, strike0, dip0, 1);
							DCFScmbopt(DCFS, temp->which_main, crst);
						}
						else okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, DCFS[temp->which_main], eqkfm2, crst, strike0, dip0, 0);
					}
					else {
						if (vary_recfault!=2) resolve_DCFS(DCFS[temp->which_main], crst, strike0, dip0, NULL, 1);
						else DCFScmbopt(DCFS, temp->which_main, crst);	//NB this does not take into account stress from aftershocks or afterslip, assuming that from mainshocks is much larger this is ok.
					}
					if (gridpoints_err==1) smoothen_DCFS(DCFS[temp->which_main], crst.nLat, crst.nLon, crst.nD, seed, 0, nn);
				}
				else {
					if (gridpoints_err==1) smoothen_DCFS(DCFS[temp->which_main], crst.nLat, crst.nLon, crst.nD, seed, 1, nn);
				}

				NFsofar+=temp->NF;
				temp=temp->next;
			}
		}
	}
}

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

void smoothen_vector(int NgridT, int nLat, int nLon, int nD, double *values, long *seed, int **nn, int return_range){
	double randcmb;
	double **interp_DCFS=dmatrix(1,NgridT, 1, 2);

	interp_nn(NgridT,nLat,nLon,nD,values,interp_DCFS,0,nn);
	for (int i=1; i<=NgridT; i++) {
		if (return_range) values[i]=interp_DCFS[i][2]-interp_DCFS[i][1];
		else {
			randcmb=interp_DCFS[i][1]+ran1(seed)*(interp_DCFS[i][2]-interp_DCFS[i][1]);
			*seed=-*seed;
			values[i]=randcmb;
		}
	}
	free_dmatrix(interp_DCFS,1,NgridT, 1, 2);
}
