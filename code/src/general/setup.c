/*
 * setup_eqkfm.c
 *
 *  Created on: Jul 24, 2013
 *      Author: camcat
 */

#include "setup.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../geom/coord_trafos.h"
#include "../geom/top_of_slipmodel.h"
#include "../inp_out/read_eqkfm.h"
#include "../inp_out/read_focmec.h"
#include "../inp_out/read_zmap.h"
#include "../okada/okadaDCFS.h"
#include "../util/error.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "../util/splines_eqkfm.h"
#include "eqkfm_copy.h"
#include "find_timesteps.h"
#include "mem_mgmt.h"
#include "struct_conversions.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int setup_catalogetc(char *catname, char **focmeccat, int nofmcat,
					 struct tm reftime, double dDCFS, double Mag_source, double Mag_main, struct crust crst,
					 struct catalog *cat, struct eqkfm **eqkfm1, double ***focmec,
					 int **firstelements, struct flags flag, int *NFM, int *Ntot,
					 double dt, double dM, double xytoll, double ztoll, double dR, double tw,
					 double tstart, double tend) {

/*	Reads catalog and focal mechanisms catalog, and fills cat, eqkfm structures accordingly.
 *
 * Input:
 *
 * catname:	ZMAP catalog file
 * focmeccat:	list of focal mechanisms catalog files [0...nofmcat-1]
 * nofmcat:	number of focal mechanisms catalogs
 * reftime: IssueTime (times will be calcualted with reference to this time)
 * dDCFS:	min. value for which grid points should be selected for calculating stress changes from source events
 * Mag_main:	magnitude of mainshocks (i.e. events which are included as sources also if flags.aftershocks==0)
 * crst:	structure containing domain information
 * flag:	flags structure
 * dt, dM, xytoll, ztoll: max. expected difference (tolerance) between events from catalog and focmec catalog;
 * dR: extra distance to be considered for sources.
 * tw:	time window to be ignored for event selection after each mainshock (NB only ignored in cat, still included as sources in eqkfm).
 * tstart: start time for including sources and catalog events //todo they should be different to include foreshocks?
 * tend: end time for forecast (sourced only included up to IssueTime, i.e. t=0)
 *
 * Output:
 *
 * cat:	catalog structure used for LL calculations
 * eqkfm1:	eqkfm structure containing stress sources. Can be NULL, and will be ignored.	[0...Ntot-1]
 * focmec:	array containing focal mechanisms [1...NFM]
 * first_elements:	indices of focmec elements which correspond to the first element of a new focal mechanism area (i.e. a new focal mechanisms catalog)
 * NFM:	length of focmec
 * Ntot: length of eqkfm1
 * Nmain:	number of mainshocks in eqkfm1
 *
 */

//todo make sure that first focal mechanism is the best oriented.

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct eqkfm *eqkfm1fm;
	double tstartS, tendS, tstartCat, tendCat;
	double minmag;
	int err=0, errP, NgridT=crst.N_allP;
	int Nfm;

	print_screen("Setting up catalog...\n");

	tendS=0;		//this is the "IssueDate", up to which data is available.
	tendCat=tend;	//this is (presumably) the "ForecastDate", up to which data is available. Events after t=0 used to calculate LL of forecast.

	//select events within some tolerance level, since they will have to be matched with focal mechanisms.
	err += readZMAP(cat, eqkfm1, Ntot, catname, crst, reftime, tstart, tendS, tstart, tendCat,
			Mag_main, tw, fmax(xytoll, dR), fmax(ztoll, dR), dDCFS, 1);

	if (err) return (err);

	// read catalog of focal mechanisms and merge it with eqkfm structure from catalog:

	if (focmeccat && (!flag.sources_all_iso || focmec)){
		err+=readmultiplefocmec(focmeccat, nofmcat, crst,fmax(xytoll, dR), fmax(ztoll, dR), dDCFS,
				reftime, tstart, tendS, tendS, (*cat).Mc, focmec, firstelements, NFM, &Nfm,  &eqkfm1fm, 1, 1);
		combine_eqkfm(*eqkfm1, eqkfm1fm, *Ntot, Nfm, dt, dM, xytoll, 1);
	}

	// filter eqkfm according to magnitude, depth.
	minmag= Mag_source;
	eqk_filter(eqkfm1, Ntot, minmag , crst.depmax+fmax(dR,ztoll));

	// calculate distances between source events and grid points.
	eqkfm2dist((*eqkfm1), crst.lat, crst.lon, crst.depth, NgridT, *Ntot, 1);


	print_screen("%d events used for catalog, %d events used as sources.\n", (int) (*cat).Z, *Ntot);
	print_logfile("%d events used for catalog, %d events used as sources.\n", (int) (*cat).Z, *Ntot);

	//todo warning if no events are selected as sources.

	return (err!=0);
}

int setup_afterslip_eqkfm(struct slipmodels_list list_slipmodels, struct crust crst, struct eqkfm **eqkfm0res){
/*
 * Reads afterslip files into eqkfm structure.
 *
 * Input:
 * 	list_slipmodel: list of slip model files
 * 	crst: structure containing crust setup information
 * 	models are the models to be used at a given time (NB: only one model per event).
 *is_afterslip indicates that all models have same geometry: Nfaults and no_slipmodels only have 1 element.
 */

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int Nm;
	double *tsnap=list_slipmodels.tsnap;
	double *tmain=list_slipmodels.tmain;
	char **slipmodels=list_slipmodels.slipmodels;
//	int *no_slipmodels=list_slipmodels.no_slipmodels;
	double *disc=list_slipmodels.disc;
	int *Nfaults=list_slipmodels.Nfaults;
    int err=0;
    int counter=0;
    int nsm=0, totfaults=0;
    char *cmb_format=list_slipmodels.cmb_format;

    //Find tot. no. of faults:
    for (int N=0; N<list_slipmodels.NSM; N++){

		Nm=list_slipmodels.no_slipmodels[N];	//number of snapshots for current afterslip event.

		if (!(strcmp(cmb_format,"farfalle"))) err+=read_farfalle_eqkfm(slipmodels[counter], NULL, Nfaults+N);
		else {
			if (!(strcmp(cmb_format,"pscmp"))) err+=read_pscmp_eqkfm(slipmodels[counter], NULL, Nfaults+N);
			else {
				if (!(strcmp(cmb_format,"fsp"))) err+=read_fsp_eqkfm(slipmodels[counter], NULL, Nfaults+N);
				else {
					print_logfile("Unknown slip model format %s (setup_afterslip_eqkfm).\n", cmb_format);
					return 1;
				}
			}
		}
		counter+=Nm;
		totfaults+=Nfaults[N];
    }

    *eqkfm0res=eqkfm_array(0,totfaults-1);
    counter=0;
    totfaults=0;
    for (int N=0; N<list_slipmodels.NSM; N++){
    	err+=setup_afterslip_element(*eqkfm0res+totfaults, slipmodels+counter, cmb_format, Nm, crst.mu, disc[N], tmain[N], tsnap+counter, crst.N_allP, crst.list_allP, list_slipmodels.cut_surf[counter], Nfaults+N, crst.lat0, crst.lon0);
		counter+=Nm;
		totfaults+=Nfaults[N];
    }

    return(err);
}


int setup_afterslip_element(struct eqkfm *eqkfm0res, char **slipmodels, char *cmb_format, int no_snap,
						double mu, double disc, double tmain, double *tsnap, int nsel,
						int *sel_pts, int cuts_surf,
						double lat0, double lon0) {

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int err=0, NF, nfmax=0;
	double 	toll=1e-10, discx, discy;
	struct eqkfm *eqkfm0;
	double ***allslip_str_temp,***allslip_dip_temp;

	struct set_of_models setmodels;

	err=read_eqkfm(slipmodels[no_snap-1], cmb_format, &eqkfm0, &NF, NULL, mu);	//find NF and eqkfm0[x].np_x
	if (err) return (err);

	//A single slip model is used for afterslip:
	setmodels.NF_models=ivector(1,1);
	setmodels.Nmod=1;
	setmodels.current_model=1;
	setmodels.NF_models[1]=NF;
	setmodels.NFmax=NF;
	setmodels.set_of_eqkfm=eqkfm_array(0, NF-1);

	// allocate temporary storage (since eqkfm0 gets overwritten)
	allslip_str_temp=malloc(NF*sizeof(double **));
	allslip_dip_temp=malloc(NF*sizeof(double **));
	for (int nf=0; nf<NF; nf++) {
		allslip_str_temp[nf]=dmatrix(0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
		allslip_dip_temp[nf]=dmatrix(0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
	}

	//read in values for slip:
	for (int m=0; m<no_snap; m++){
		err=read_eqkfm(slipmodels[m], cmb_format, &eqkfm0, &NF, NULL, mu);
		for (int nf=0; nf<NF; nf++){
			eqkfm0[nf].tot_slip[m]=eqkfm0[nf].tot_slip[0];
			copy_vector(eqkfm0[nf].slip_str, &(allslip_str_temp[nf][m]), eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			copy_vector(eqkfm0[nf].slip_dip, &(allslip_dip_temp[nf][m]), eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			free_dvector(eqkfm0[nf].slip_str,1,0);
			free_dvector(eqkfm0[nf].slip_dip,1,0);
		}

		if (err) return (err);
	}

	if (cuts_surf) {
			top_of_slipmodel(eqkfm0, NF);
			for (int i=0; i<NF; i++) eqkfm0[i].cuts_surf=1;
	}

	for (int nf=0; nf<NF; nf++) {
		eqkfm0[nf].tot_slip=dvector(0,no_snap-1);
		eqkfm0[nf].ts=dvector(1,no_snap);	//shifted by one element because of indexing in copy_vector function.
		copy_vector(tsnap-1, &(eqkfm0[nf].ts), no_snap);	//copy afterslip time steps into eqkfm0 structure.
		eqkfm0[nf].ts+=1;	//since should start from 0th element (not 1st).
		eqkfm0[nf].nosnap=no_snap;
		eqkfm0[nf].t=tmain;
		eqkfm0[nf].nsel=nsel;
		eqkfm0[nf].selpoints=sel_pts;
		//todo check: this line ok? (copied from CRSjuly)
		eqkfm0[nf].is_slipmodel=1;
		latlon2localcartesian(eqkfm0[nf].lat, eqkfm0[nf].lon, lat0, lon0, &(eqkfm0[nf].y), &(eqkfm0[nf].x));
		eqkfm0[nf].allslip_str=allslip_str_temp[nf];
		eqkfm0[nf].allslip_dip=allslip_dip_temp[nf];
		setmodels.set_of_eqkfm[nf]=eqkfm0[nf];
	}

	//allocate memory and copy values from setmodels;
	(*eqkfm0res).parent_set_of_models=(struct set_of_models *) malloc((size_t) (sizeof(struct set_of_models)));
	memcpy((*eqkfm0res).parent_set_of_models, &setmodels, (size_t) sizeof(struct set_of_models));

	set_current_slip_model(eqkfm0res,1);
	
	return err;
}



int setup_eqkfm_element(struct eqkfm *eqkfm0res, char **slipmodels, char *cmb_format, int no_slipmodels,
						double mu, double disc, double tmain, int nsel,
						int *sel_pts, double *mmain, int cuts_surf,
						int *NF0, double lat0, double lon0) {

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct set_of_models setmodels;
	struct eqkfm *eqkfm0;
	int err=0, NF, nftot=0, nfmax=0;
	double 	toll=1e-10, discx, discy;

	setmodels.NF_models=ivector(1,no_slipmodels);
	setmodels.Nmod=no_slipmodels;
	setmodels.current_model=1;

	for (int m=0; m<no_slipmodels; m++){
		err=read_eqkfm(slipmodels[m], cmb_format, NULL, &NF, NULL, mu);
		if (err){
			print_screen(" ** Error: Input slip model %s could not be read (setup_eqkfm_element). **\n", slipmodels[m]);
			print_logfile(" ** Error: Input slip model %s could not be read (setup_eqkfm_element). **\n", slipmodels[m]);
			return (1);
		}
		setmodels.NF_models[m+1]=NF;
		nfmax=fmax(nfmax,NF);
		nftot+=NF;
	}
	setmodels.NFmax=nfmax;
	setmodels.set_of_eqkfm=eqkfm_array(0,nftot-1);

	nftot=0;
	for (int m=0; m<no_slipmodels; m++){
		err=read_eqkfm(slipmodels[m], cmb_format, &eqkfm0, &NF, mmain, mu);
		if (err) return (err);

		if (cuts_surf) {
			top_of_slipmodel(eqkfm0, NF);
			for (int i=0; i<NF; i++) eqkfm0[i].cuts_surf=1;
		}

		for (int nf=0; nf<NF; nf++) {
			eqkfm0[nf].t=tmain;
			eqkfm0[nf].nsel=nsel;
			eqkfm0[nf].selpoints=sel_pts;
			//todo check: this line ok? (copied from CRSjuly)
			eqkfm0[nf].is_slipmodel=1;
			latlon2localcartesian(eqkfm0[nf].lat, eqkfm0[nf].lon, lat0, lon0, &(eqkfm0[nf].y), &(eqkfm0[nf].x));
			setmodels.set_of_eqkfm[nftot+nf]=eqkfm0[nf];
		}
		nftot+=NF;
	}

	//allocate memory and copy values from setmodels;
	(*eqkfm0res).parent_set_of_models=(struct set_of_models *) malloc((size_t) (sizeof(struct set_of_models)));
	memcpy((*eqkfm0res).parent_set_of_models, &setmodels, (size_t) sizeof(struct set_of_models));

	set_current_slip_model(eqkfm0res,1);
	if (NF0) *NF0=nfmax;
	return err;
}

void set_current_slip_model(struct eqkfm *eqkfm0, int slipmodel_index){
//sets variables in eqkfm0 to required slip model.

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int nftot=0;
	struct set_of_models allmod= *(eqkfm0[0].parent_set_of_models);
	struct eqkfm* eqkfmall=allmod.set_of_eqkfm;

	// allmod.Nmod==0 means that no slip model is used for this event: do nothing.
	if (!allmod.Nmod) return;

	// find where slip model no.slipmodel_index starts;
	// copy slip models from eqkfmall to eqkfm0;
	// delete all info from remaining elements of eqkfm0;
	for (int i=1; i<slipmodel_index; i++) nftot+=allmod.NF_models[i];
	for (int nf=0; nf<allmod.NF_models[slipmodel_index]; nf++) copy_eqkfm_all(eqkfmall[nf+nftot], eqkfm0+nf);
	for (int nf=allmod.NF_models[slipmodel_index]; nf<allmod.NFmax; nf++) empty_eqkfm(eqkfm0+nf);

	print_logfile("Slip model set to no. %d.\n", slipmodel_index);

	return;
}

int setup_CoeffsDCFS(struct Coeff_LinkList **Coefficients, struct pscmp **DCFS_out,
		struct crust crst, struct eqkfm *eqkfm0, int Nm, int *Nfaults, double aftersliptime, int afterslip) {
	/*
	 *  aftersliptime: event time of the mainshock containing afterslip.
	 *  int afterslip: flag indicating if afterslip should be used.
	 */

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct pscmp *DCFS;
    struct Coeff_LinkList *AllCoeff, *temp;
    int NFsofar=0, Nsel, Nsteps, NFtot, eq;
    int mainshock_withafterslip;
    double M0;

    //----------set up Coefficients----------------//

    if (Coefficients!=NULL) {
    	//todo make sure that coefficients are only recalculated when needed (i.e. only for events for which several slip model are available).

    	//Create elements of structure (allocate memory):
    	AllCoeff= malloc( sizeof(struct Coeff_LinkList));	//TODO deallocate at the end.
		temp= AllCoeff;
		for(int i=0; i<Nm; i++) {
			if (i<Nm-1) {
				temp->hasafterslip=0;
				temp->next= malloc(sizeof(struct Coeff_LinkList));
				temp= temp->next;
			}
			else {
				temp->hasafterslip=0;
				temp->next=(struct Coeff_LinkList *) 0;
			}
		}

		*Coefficients=AllCoeff;
		print_logfile("Okada Coefficients structure set up.\n");
    }

    //--------------set up DCFS-------------------//

    if (DCFS_out!=NULL){

		DCFS=pscmp_arrayinit(crst,0,Nm-1);

		NFtot=0;
		for (int eq=0; eq<Nm; eq++){
			Nsel = eqkfm0[NFtot].nsel;
			DCFS[eq].fdist=eqkfm0[NFtot].distance;
			DCFS[eq].index_cat=eqkfm0[NFtot].index_cat;
			DCFS[eq].which_pts=eqkfm0[NFtot].selpoints;
			DCFS[eq].t=eqkfm0[NFtot].t;
			M0=0.0;
			for (int f=0; f<Nfaults[eq]; f++) M0+=pow(10,1.5*(eqkfm0[NFtot+f].mag+6));
			DCFS[eq].m=(2.0/3.0)*log10(M0)-6;
			DCFS[eq].NF=Nfaults[eq];
			if (DCFS[eq].nsel!=Nsel){
				// todo [coverage] this block is never tested
				if (DCFS[eq].nsel>0){
					free_d3tensor(DCFS[eq].S,1,DCFS[eq].nsel,1,3,1,3);
					free_dvector(DCFS[eq].cmb,1,DCFS[eq].nsel);
				}
				DCFS[eq].nsel=Nsel;
				DCFS[eq].S = d3tensor(1,Nsel,1,3,1,3);
				DCFS[eq].cmb=dvector(1,Nsel);
				for (int i=1; i<=Nsel; i++) DCFS[eq].cmb[i]=0.0;
			}
			NFtot+=Nfaults[eq];
		}

		*DCFS_out=DCFS;
		print_logfile("DCFS structure set up.\n");
    }

    //--------------associates afterslip with one mainshock-------------------//
	// uses a ~1 sec tolerance

    if (afterslip){
    	int i=0;
		mainshock_withafterslip=closest_element(timesfrompscmp(DCFS, Nm), Nm, aftersliptime, 0.000011575);
		if (mainshock_withafterslip==-1){
			print_logfile("Error: Reference time for afterslip does not correspond to a mainshock. Exiting.\n");
			print_screen("Error: Reference time for afterslip does not correspond to a mainshock. Exiting.\n");
			return(1);
		}
		struct Coeff_LinkList *temp;
		temp=AllCoeff;
		while (i<Nm && i!=mainshock_withafterslip) {
			i++;
			temp=temp->next;
		}
		temp->hasafterslip=1;
    }

    return(0);
}

int update_CoeffsDCFS(struct Coeff_LinkList **Coefficients,
		struct crust crst, struct eqkfm *eqkfm0, int Nm, int *Nfaults) {

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

    struct Coeff_LinkList *temp;
    struct set_of_models tmp_setofmodels;
    int NFsofar=0;
    static int first_timein=1, switch_slipmodel;

	//Fill in elements of structure:
	temp= *Coefficients;
	for(int i=0; i<Nm; i++) {
		if (eqkfm0[NFsofar].is_slipmodel) {
			// coefficients should only calculated if more than one slip model is provided (switch_slipmodel):
			tmp_setofmodels=*(eqkfm0[NFsofar].parent_set_of_models);
			switch_slipmodel= (tmp_setofmodels.Nmod > 1);
			if (first_timein || switch_slipmodel){
				if (!first_timein){
					if (temp->Coeffs_st) free_f3tensor(temp->Coeffs_st, 1,0,1,0,1,0);
					if (temp->Coeffs_dip) free_f3tensor(temp->Coeffs_dip, 1,0,1,0,1,0);
				}
				okadaCoeff(&(temp->Coeffs_st), &(temp->Coeffs_dip), eqkfm0+NFsofar,
					   Nfaults[i], crst, crst.lat, crst.lon, crst.depth);
//				okadaCoeff_mpi(&(temp->Coeffs_st), &(temp->Coeffs_dip), eqkfm0+NFsofar,
//					   Nfaults[i], crst, crst.lat, crst.lon, crst.depth);
				temp->NgridT=eqkfm0[0].nsel;
				temp->NF=Nfaults[i];
				temp->NP=0;
				for(int f=0; f<Nfaults[i]; f++) {
					temp->NP+= eqkfm0[NFsofar+f].np_di*eqkfm0[NFsofar+f].np_st;
				}
				temp->which_main=i;
			}
		}
		else {
			temp->Coeffs_st=temp->Coeffs_dip=NULL;
			if (first_timein){
				temp->NgridT=eqkfm0[0].nsel;
				temp->NF=Nfaults[i];
				temp->NP=0;
				for(int f=0; f<Nfaults[i]; f++) {
					temp->NP+= eqkfm0[NFsofar+f].np_di*eqkfm0[NFsofar+f].np_st;
				}
				temp->which_main=i;
			}
		}
		NFsofar+=Nfaults[i];
		if (i<Nm-1) {
			temp= temp->next;
		}
	}
	print_logfile("Okada Coefficients structure updated.\n");

	first_timein=0;

    return(0);
}

int setup_afterslip_evol(double t0, double t1, double *Cs, double *ts,
						 int Nfun, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2,
						 long *seed) {

//if splines are, eqk_aft is substituted with more densily discretized version (L steps instead of Nas).
//Cs, ts, coefficients of temporal evolution functions (See Savage Parkfield paper). Nfun: no. of such funtions.
// eq_aft has indices: [0...Nas*Nfaults-1].

	int err=0;
	int L0=0, i;
	int splines;
	double TAU=200000, dtau=7000;	//todo allow to set from outside.
	double Kotau;
	double M0,mu;
	double smallstepstime=12;
	double now, prev, norm, curr;
	double Tendaft;	//Time to which cumulative afterslip snapshot refers.
	struct eqkfm *eq_aft;

	double Teq;
	double *t_afterslip;
	int Nas=(*eqk_aft)[0].nosnap;
	int Ltot=0, offset=1, nfaults;

	eq_aft= *eqk_aft;

	splines=(Nas>1)? 1 : 0;	//fixme for each nev
	//todo find L first, then setup vectors.
	*times2=dvector(0,*L-1);

	int nev=0;
	double tend;

	offset= (t0<(*eqk_aft)[0].t) ? 1 : 0;	//need an extra time step at the start if t0<Teq[0].
	if (offset) (*times2)[0]=t0;

	//fixme the warning in findtimestepsomori will not work (since L is overwritten, and L0 not even initialized).
	nev=nfaults=0;
	while(nev<NA && (*eqk_aft)[nfaults].t<t1){
		Teq=(*eqk_aft)[nfaults].t;
		tend= (nev<NA-1) ? fmin((*eqk_aft)[nfaults+Nfaults[nev]].t,t1) : t1;

		err+=findtimestepsomori(Teq, Teq, fmin(smallstepstime+Teq,tend), 0, 183, TAU, 0.3*dtau, 0.6, 0.001, (*times2)+offset, &Kotau, L);
		if (smallstepstime+Teq<tend) err+=findtimestepsomori(Teq, smallstepstime+Teq, tend, 0, 183, TAU, dtau, 0.6, 0.001, (*times2)+*L+offset, &Kotau, &L0);
		Ltot=offset=*L+L0+offset;
		nev++;
		nfaults+=Nfaults[nev];
	}

	*L=Ltot;

	// Temporal evolution of afterslip.//
	if (splines==0){
		nfaults=0;
		for (nev=0; nev<NA; nev++){
			Tendaft=(*eqk_aft)[nfaults].ts[(*eqk_aft)[nfaults].nosnap-1];
			Teq=(*eqk_aft)[nfaults].t;

			//allocate tevol vectors:
			for (int f=0; f<Nfaults[nev]; f++) {
				(*eqk_aft)[nfaults+f].tevol=dvector(0,*L-1);
			}

			norm=0.0;
			for (int i=0; i<Nfun; i++) norm+=Cs[i]*log(1+(Tendaft-Teq)/ts[i]);
			now=0.0;
			curr=0.0;

			for (int t=1; t<=*L; t++){
				if ((*times2)[t]<Teq){
					for (int f=0; f<Nfaults[nev]; f++) {
						(*eqk_aft)[nfaults+f].tevol[t-1]= 0.0;
					}
				}
				else{
					prev=curr;
					now= (*times2)[t]-Teq;
					curr=0.0;
					for (int i=0; i<Nfun; i++) curr+=Cs[i]*log(1+now/ts[i]);
					for (int f=0; f<Nfaults[nev]; f++) {
						(*eqk_aft)[nfaults+f].tevol[t-1]= (curr-prev)/norm;
					}
				}
			}
			nfaults+=Nfaults[nev];
		}
	}
	else{
		nfaults=0;
		for (nev=0; nev<NA; nev++){
			Teq=(*eqk_aft)[nfaults].t;
			t_afterslip=(*eqk_aft)[nfaults].ts;
			Nas=(*eqk_aft)[nfaults].nosnap;
			for (int i=0; i<Nas; i++) t_afterslip[i]-=Teq;	//since functions below start from t=0;
			for (int i=0; i<=*L; i++) (*times2)[i]-=Teq;	//since functions below start from t=0;

			splines_eqkfm(&eq_aft, Nas, Nfaults[nev], t_afterslip-1, (*times2)-1, *L, seed);

			for (int f=0; f<Nfaults[nev]; f++) {

				for (int p=1; p<=eq_aft[f].np_di*eq_aft[f].np_st; p++) {
					for (int l=*L-1; l>=0; l--) {
						if ((*times2)[l]<Teq){
							//no afterslip before its start time:
							eq_aft[f].allslip_str[l][p]=0.0;
							eq_aft[f].allslip_dip[l][p]=0.0;
						}
						else{
							if (l>0 && (*times2)[l-1]>Teq){	//this is to avoid subtracting from element with t<Teq.
								eq_aft[f].allslip_str[l][p]-=eq_aft[f].allslip_str[l-1][p];
								eq_aft[f].allslip_dip[l][p]-=eq_aft[f].allslip_dip[l-1][p];
							}
						}
					}
				 }
				eq_aft[f].tevol=NULL;
			}
			for (int i=0; i<Nas; i++) t_afterslip[i]+=Teq;	//revert shift from before;
			for (int i=0; i<=*L; i++) (*times2)[i]+=Teq;	//revert shift from before;
			nfaults+=Nfaults[nev];
			eq_aft+=Nfaults[nev];	//gets shifted every time so 0th element is the one of the new nev.
		}

		*eqk_aft=eq_aft-nfaults;	//shift it back.
	}

	return(err!=0);
}
