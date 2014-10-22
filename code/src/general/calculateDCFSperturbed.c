/*
 * calculateDCFSperturbed.c
 *
 *  Created on: Mar 19, 2013
 *      Author: camcat
 */

#include "calculateDCFSperturbed.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

void calculateDCFSperturbed(double **DCFSrand, struct pscmp *DCFS, struct eqkfm *eqkfmAf,
							struct eqkfm *eqkfm0, struct flags flag,
							double *tevol, double *times, int Nmain, struct crust crst,
							struct Coeff_LinkList *AllCoeff, int NTScont,
							double **focmec, int *fmzoneslim, int NFM, long *seed,
							double tdata0, double tdata1, int refresh, int which_recfault) {

/* Input:
 * 
 * eqkfmAf contains afterslip snapshots;
 * eqkfm0 contains seismic sources;
 * tevol: if continuous process is stationary (splines==0), tevol contains time dependence of the process;
 * times: times to which elements of tevol correspond: S(t=times[j])=S0*tevol[j], where S=slip, S0 is the slip contained in eqkfmAf.
 * NTScont is the total number of time steps for continuous process;
 * Nmain: length of eqkfm0;	//todo change variable name
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
 * NB: not general, as it assumes mainshock is the first event. todo check/fix this.
 */

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	//flags:
	int	afterslip=flag.afterslip, \
		vary_recfault=flag.err_recfault, \
		new_slipmodel=flag.new_slipmodel, \
		gridpoints_err=flag.err_gridpoints, \
		splines=flag.splines, \
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
		print_logfile("\nSetting up variables for calculating perturbed Coulomb fields.\n");

		new_slipmodel=1;
		//-----calculate neighbouring points--------//

		if (gridpoints_err==1) {
			nn=imatrix(1,NgridT,1,6);
			nearest_neighbours(NgridT, crst.nLat,crst.nLon,crst.nD, nn);
			interp_DCFS=dmatrix(1,NgridT,1,2);
			mycmb=dvector(1,NgridT);
		}

		if (crst.variable_fixmec || (vary_recfault==1 && which_recfault==0)){
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

			Coeffs_st=AllCoeff->Coeffs_st;
			Coeffs_dip=AllCoeff->Coeffs_dip;
			for (int i=0; i<DCFS_Af_size; i++)	{
				okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, DCFS_Af[i], eqkfmAf+i*DCFS_Af[0].NF, crst, NULL, NULL, 1); //todo free memory used by *AllCoeff; todo make this work for splines==1 too...
				if (vary_recfault==0) {
					resolve_DCFS(DCFS_Af[i], crst, crst.str0, crst.dip0, NULL, 1);
					if (splines && i<DCFS_Af_size-1){
						for (int n=1; n<=NgridT; n++) cmb_cumu[n]+=DCFS_Af[i].cmb[n];
					}
				}
			}
		}
	}


		//-----------------------------------------------------------------//
		//							Mainshock							   //
		//-----------------------------------------------------------------//

	if (time_in==1 || refresh){
		if (afterslip!=2){
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
		print_screen("*Error: function calculateDCFSperturbed doesn't know how to calculate OOPS when afterslip is included!!*\n");
		print_logfile("*Error: function calculateDCFSperturbed doesn't know how to calculate OOPS when afterslip is included!!*\n");
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
			}
			else {
		// todo [coverage] this block is never tested
				*strike0=focmec[1][which_recfault];
				*dip0=focmec[2][which_recfault];
				*rake0=focmec[3][which_recfault];	//only used for splines==1 (see below).
			}
			break;

		case 0:
			if (crst.variable_fixmec){
				strike0= crst.str0;
				dip0=crst.dip0;
			}

			else{
				*strike0= crst.str0[0];
				*dip0=crst.dip0[0];
			}
			*rake0=crst.rake0[0];	//only used for splines==1 (see below). todo allow rake to be passed by user as well.
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
		// todo [coverage] this line is never tested
			for (int n=1; n<=NgridT; n++) if (DCFSrand) DCFSrand[l][n]=0.0;	//todo should DCFSrand just set to NULL if afterslip==0?
		}
	}
	else {
		if (times[0]<tdata1 && times[NTScont-1]>=tdata0){

			if (splines==0){
				if (vary_recfault!=0) resolve_DCFS(DCFS_Af[0], crst, strike0, dip0, NULL, 1);
				for (int l=0; l<NTScont; l++) {
					if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
					for (int n=1; n<=NgridT; n++) DCFSrand[l][n]=tevol[l]*DCFS_Af[0].cmb[n];
				}
			}

			else{
				if (vary_recfault!=0){
					// todo [coverage] this block is never tested
					for (int n=1; n<=NgridT; n++) cmb_cumu[n]=0.0;
					for (int l=0; l<NTScont; l++) {
						resolve_DCFS(DCFS_Af[l], crst, strike0, dip0, NULL, 1);
						if (l<NTScont-1) for (int n=1; n<=NgridT; n++) cmb_cumu[n]+=DCFS_Af[l].cmb[n];
					}
				}

				for (int l=0; l<NTScont; l++) {
					if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
					for (int n=1; n<=NgridT; n++) {
						DCFSrand[l][n]= (fabs(cmb_cumu[n])>DCFS_cap) ? (DCFS_cap/fabs(cmb_cumu[n]))*DCFS_Af[l].cmb[n] : DCFS_Af[l].cmb[n];
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
					Coeffs_st=temp->Coeffs_st;
					Coeffs_dip=temp->Coeffs_dip;
					if (vary_recfault!=2) resolve_DCFS(DCFS[temp->which_main], crst, strike0, dip0, NULL, 1);
					// todo [coverage] this block is never tested	
					else DCFScmbopt(DCFS, temp->which_main, crst);	//NB this does not take into account stress from afterslip, assuming that from mainshocks is much larger this is ok.

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
