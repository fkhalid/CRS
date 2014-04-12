/*
 * setup_eqkfm.c
 *
 *  Created on: Jul 24, 2013
 *      Author: camcat
 */

#include "setup.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int setup_catalogetc(char *catname, char **focmeccat, int nofmcat, char *fm_format,
					 struct tm reftime, double dDCFS, double Mag_main, struct crust crst,
					 struct catalog *cat, struct eqkfm **eqkfm1, double ***focmec,
					 int **firstelements, struct flags flag, int *NFM, int *Ntot, int *Nmain,
					 double dt, double dM, double xytoll, double ztoll, double dR, double tw,
					 double tstart, double tend, double Mc_source) {
//  focmec will contain matrix with focal mechanisms (*NFM of them). can be null, and will be ignored.
//	double t0, double tw1, double tw2, double t1 =	tstartLL, tmain[0], tmain[0]+tw, fmax(tendLL,Tend)
//  dt, dM, xytoll, ztoll: expected difference b/t/ events from catalog and focmec catalog; dR= extra distance to be considered for sources.
//	focal mechanisms are ignored if focmeccat is NULL.
//  tstartLL: start of LL inversion (only focal mechanisms before this time will be included).
//  tstartData: time after which which events should be included as sources.


//todo make sure that first focal mechanism is the best oriented.
//todo only include as sources events large enough to affect more than 1 grid point...
//todo so far, many events are selected and  eqkfm[n].nsel=0 is used to deactivate them. This is a waste (done this way to avoid messing up indices, but should be changed).
//todo check Mc_source works.

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct eqkfm *eqkfm1fm;
	double tstartS, tendS, tstartCat, tendCat;
	int err=0, errP, NgridT=crst.N_allP;
	int Nfm;
	int numCols;

	if(procId == 0) {
		if (verbose_level>0) printf("Setting up catalog...\n");
	}

	tendS=0;		//this is the "IssueDate", up to which data is available.
	tendCat=tend;	//this is (presumably) the "ForecastDate", up to which data is available. Events after t=0 used to calculate LL of forecast.

	//select events within some tolerance level, since they will have to be matched with focal mechanisms.
	//todo remove this line (only for making it work for Parkfield).

	if(procId == 0) {
		numCols = countcol(catname);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&numCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(numCols == 8) {
		err += read_RS(catname, cat, crst, -2.0, tstart, 0.0, tw, tendCat, 0.0, eqkfm1, dDCFS, Ntot, 1);
		for (int eq=0; eq<*Ntot; eq++) {
			if ((*eqkfm1)[eq].t>0.0) {
				*Ntot=eq;
				break;
			}
		}
	}
	else {
		err += readZMAP(cat, eqkfm1, Ntot, catname, crst, reftime, tstart, tendS, tstart, tendCat,
						Mag_main, tw, fmax(xytoll, dR), fmax(ztoll, dR), dDCFS, 1);
	}

	//fixme check that foc mec are read when aftershocks==1.
	//if (flag.aftershocks){
		//if (focmec){
			err+=readmultiplefocmec(focmeccat, nofmcat, fm_format, crst,fmax(xytoll, dR), fmax(ztoll, dR), dDCFS,
					reftime, tstart, tendS, tendS, (*cat).Mc, focmec, firstelements, NFM, &Nfm,  &eqkfm1fm, 1, 1);
			//err+=readfocmec(focmeccat[0], fm_format, crst, fmax(xytoll, dR), fmax(ztoll, dR), dDCFS, reftime, tstartS, tendS, tendS, (*cat).Mc, focmec, NFM, &Nfm, &eqkfm1fm, 1, 1);
			errP=combine_eqkfm(*eqkfm1, eqkfm1fm, *Ntot, Nfm, tendS, dt, dM, xytoll, 1);
			//err+=(errP==NULL);	//commented since foc mec catalog may not be available.
		//}
		if (flag.aftershocks) eqk_filter(eqkfm1, Ntot, (Mc_source>20) ? (*cat).Mc : Mc_source, crst.depmax+fmax(dR,ztoll));
		else eqk_filter(eqkfm1, Ntot, Mag_main, crst.depmax+fmax(dR,ztoll));	//only keep mainshocks.

		eqkfm2dist((*eqkfm1), crst.lat, crst.lon, crst.depth, NgridT, *Ntot, 1);
	//}

//	else {
//		if (focmec){
//			err+=readmultiplefocmec(focmeccat, nofmcat, fm_format, crst,fmax(xytoll, dR), fmax(ztoll, dR), dDCFS,
//					reftime, tstart, tendS, tendS, Mag_main, focmec, firstelements, NFM, &Nfm,  &eqkfm1fm, 1, 1);
//			errP=combine_eqkfm(*eqkfm1, eqkfm1fm, *Ntot, Nfm, tendS, dt, dM, xytoll, 1);
//		}
//		eqk_filter(eqkfm1, Ntot, Mag_main, crst.depmax+fmax(dR,ztoll));	//only keep mainshocks.
//		eqkfm2dist((*eqkfm1), crst.lat, crst.lon, crst.depth, NgridT, *Ntot, 0);
//	}

	if (Nmain) *Nmain=0;
	for (int i=0; i<(*Ntot); i++) {
		if ((*eqkfm1)[i].mag>=Mag_main) {
			(*eqkfm1)[i].is_mainshock=1;
			if (Nmain) *Nmain+=1;
		}
		else{
			if (flag.only_aftershocks_withfm && !(*eqkfm1)[i].is_slipmodel) (*eqkfm1)[i].nsel=0;
			if (flag.full_field==0) (*eqkfm1)[i].is_slipmodel=0;
		}
	}

	if(procId == 0) {
		if (verbose_level>0) printf("%d events used for catalog, %d events used as sources, %d of which mainshocks.\n", (int) (*cat).Z, *Ntot, *Nmain);
		if (flog) {
			fprintf(flog,"%d events used for catalog, %d events used as sources, %d of which mainshocks.\n", (int) (*cat).Z, *Ntot, *Nmain);
			fflush(flog);
		}
	}

	return (err!=0);
}

int setup_afterslip_eqkfm(struct slipmodels_list list_slipmodels, struct crust crst, int resample, struct eqkfm **eqkfm0res){
//input models are the models to be used at a given time (NB: only one model per event).
//is_afterslip indicates that all models have same geometry: Nfaults and no_slipmodels only have 1 element.

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int Nm=list_slipmodels.NSM;
	double *tmain=list_slipmodels.tmain;
	char **slipmodels=list_slipmodels.slipmodels;
	int *no_slipmodels=list_slipmodels.no_slipmodels;
	double *disc=list_slipmodels.disc;
	int *Nfaults=list_slipmodels.Nfaults;
	double d_touching_faults=3.0;	//todo: set somewhere else.
    int NFtot=0;
    int err=0;


    	if (!(strcmp(cmb_format,"farfalle"))) err+=read_farfalle_eqkfm(slipmodels[0], NULL, Nfaults);
    	else {
			if (!(strcmp(cmb_format,"pscmp"))) err+=read_pscmp_eqkfm(slipmodels[0], NULL, Nfaults);
			else {
				if (!(strcmp(cmb_format,"fsp"))) err+=read_fsp_eqkfm(slipmodels[0], NULL, Nfaults);
				else {
					if(procId == 0) {
						if (flog) {
							fprintf(flog,"Unknown slip model format %s (setup_afterslip_eqkfm).\n", cmb_format);
							fflush(flog);
						}
					}

					return 1;
				}
			}
    	}
    	NFtot=Nfaults[0]*Nm;

    *eqkfm0res=eqkfm_array(0, NFtot-1);

    //TODO: implement multiple slip models with non strike slip event.
    NFtot=0;
    for (int nn=0; nn<Nm; nn++){
    	err+=setup_eqkfm_element((*eqkfm0res)+NFtot, slipmodels+nn, no_slipmodels[0], crst.mu, disc[0], tmain[nn], d_touching_faults, crst.N_allP, crst.list_allP, NULL, resample, 0, Nfaults, crst.lat0, crst.lon0);
		NFtot+=Nfaults[0];
	}

    return(err);
}

int setup_eqkfm_element(struct eqkfm *eqkfm0res, char **slipmodels, int no_slipmodels,
						double mu, double disc, double tmain, double d_close, int nsel,
						int *sel_pts, double *mmain, int resample, int tap_bot, int *NF0,
						double lat0, double lon0) {
	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	static struct set_of_models setmodels;
	static struct eqkfm *eqkfmall;
	struct eqkfm *eqkfm0;
	int err=0, NF, nftot=0, nfmax=0;
	double 	toll=1e-10, discx, discy;

	(*eqkfm0res).parent_set_of_models=&setmodels;
	setmodels.NF_models=ivector(1,no_slipmodels);
	setmodels.Nmod=no_slipmodels;
	setmodels.current_model=1;

	for (int m=0; m<no_slipmodels; m++){
		err=read_eqkfm(slipmodels[m], NULL, &NF, NULL, mu);
		if (err){
			if(procId == 0) {
				if (verbose_level>0) printf(" ** Error: Input slip model %s could not be read (setup_eqkfm_element). **\n", slipmodels[m]);
				if (flog) fprintf(flog, " ** Error: Input slip model %s could not be read (setup_eqkfm_element). **\n", slipmodels[m]);
			}

			return (1);
		}
		setmodels.NF_models[m+1]=NF;
		nfmax=fmax(nfmax,NF);
		nftot+=NF;
	}
	setmodels.NFmax=nfmax;
	eqkfmall=eqkfm_array(0,nftot-1);
	setmodels.set_of_eqkfm=eqkfmall;

	nftot=0;
	for (int m=0; m<no_slipmodels; m++){
		err=read_eqkfm(slipmodels[m], &eqkfm0, &NF, mmain, mu);
		if (err) return (err);
		which_taper(eqkfm0, NF, tap_bot, 0, d_close);
		for (int nf=0; nf<NF; nf++) {
			eqkfm0[nf].t=tmain;
			eqkfm0[nf].nsel=nsel;
			eqkfm0[nf].selpoints=sel_pts;
			latlon2localcartesian(eqkfm0[nf].lat, eqkfm0[nf].lon, lat0, lon0, &(eqkfm0[nf].y), &(eqkfm0[nf].x));
			if (resample) {
			  	discx=eqkfm0[nf].L/eqkfm0[nf].np_st;
			  	discy=eqkfm0[nf].W/eqkfm0[nf].np_di;
			  	//if current discretization is larger than required, resample:
			  	if (discx>disc || discy>disc) {
			  		suomod1_resample(eqkfm0[nf], eqkfmall+nftot+nf, disc, 0.0);
					if(procId == 0) {
				  		if (flog){
							fprintf(flog, "slip model %s is resampled from res=[str=%.3lf, dip=%.3lf] to res=%.3lf (setup.c).\n", slipmodels[m], discx, discy, disc);
							fflush(flog);
						}
					}
			  	}
			  	else {
					if (fabs(discx-discy)>toll) {
						err+=suomod1_resample(eqkfm0[nf], eqkfmall+nftot+nf, fmin(discx, discy), 0.0);	//create a slip model with square patches.
						if (flog){
							fprintf(flog, "slip model %s is resampled to obtain square patches (setup.c).\n", slipmodels[m]);
							fflush(flog);
						}
					}
					else copy_eqkfm_all(eqkfm0[nf], eqkfmall+nftot+nf);
			  	}
				err+=suomod1_taper(eqkfmall[nftot+nf], eqkfmall+nftot+nf);	//just taper old slip model.
			}
			else err+=suomod1_taper(eqkfm0[nf], eqkfmall+nftot+nf);	//just taper old slip model.
		}
		nftot+=NF;
	}

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

	if(procId == 0) {
		if (flog) {
			fprintf(flog,"Slip model set to no. %d.\n", slipmodel_index);
			fflush(flog);
		}
	}

	return;
}

int setup_CoeffsDCFS(struct Coeff_LinkList **Coefficients, struct pscmp **DCFS_out,
		struct crust crst, struct eqkfm *eqkfm0, struct eqkfm *eqkfm1, int Nm,
		int Ntot, int *Nfaults, int *which_main) {
	//set Ntot=0 if aftershocks should not be considered.

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct pscmp *DCFS;
    struct Coeff_LinkList *AllCoeff, *temp;
    int NFsofar=0, Nsel, Nsteps, NFtot, eq;
    double M0;

    //----------set up Coefficients----------------//

    if (Coefficients!=NULL) {
		AllCoeff= malloc( sizeof(struct Coeff_LinkList));	//TODO deallocate at the end.
		temp= AllCoeff;
		for(int i=0; i<Nm; i++) {
			if (eqkfm0[NFsofar].is_slipmodel) {
				okadaCoeff(&(temp->Coeffs_st), &(temp->Coeffs_dip), eqkfm0+NFsofar,
						   Nfaults[i], crst, crst.lat, crst.lon, crst.depth);
//				okadaCoeff_mpi(&(temp->Coeffs_st), &(temp->Coeffs_dip), eqkfm0+NFsofar,
//						   Nfaults[i], crst, crst.lat, crst.lon, crst.depth);
			}
			else {
				temp->Coeffs_st=temp->Coeffs_dip=NULL;
			}
			temp->NgridT=eqkfm0[0].nsel;
			temp->NF=Nfaults[i];
			temp->NP=0;
			for(int f=0; f<Nfaults[i]; f++) {
				temp->NP+= eqkfm0[NFsofar+f].np_di*eqkfm0[NFsofar+f].np_st;
			}
			temp->which_main=which_main[i];
			NFsofar+=Nfaults[i];
			if (i<Nm-1) {
				temp->next= malloc(sizeof(struct Coeff_LinkList));
				temp= temp->next;
			}
			else temp->next=(struct Coeff_LinkList *) 0;
		}
		*Coefficients=AllCoeff;
		if(procId == 0) {
			if (flog){
				fprintf(flog,"Okada Coefficients structure set up.\n");
				fflush(flog);
			}
		}
    }

    //--------------set up DCFS-------------------//

    if (DCFS_out!=NULL){
		Nsteps= fmax(Ntot,Nm);
		DCFS=pscmp_arrayinit(crst,0,Nsteps-1);

		for (int eq1=0; eq1<Ntot; eq1++){
			Nsel = eqkfm1[eq1].nsel;
			DCFS[eq1].NF=1;
			DCFS[eq1].index_cat=eqkfm1[eq1].index_cat;
			DCFS[eq1].which_pts=eqkfm1[eq1].selpoints;
			DCFS[eq1].fdist=eqkfm1[eq1].distance;
			DCFS[eq1].nsel=Nsel;
			DCFS[eq1].t=eqkfm1[eq1].t;
			DCFS[eq1].m=eqkfm1[eq1].mag;
			if (Nsel>0){
				DCFS[eq1].S = d3tensor(1,Nsel,1,3,1,3);
				DCFS[eq1].cmb=dvector(1,Nsel);
				for (int i=1; i<=Nsel; i++) DCFS[eq1].cmb[i]=0.0;
			}
		}

		NFtot=0;
		//substitute mainshocks:
		for (int eq1=0; eq1<Nm; eq1++){
			eq= which_main[eq1];
			Nsel = eqkfm0[NFtot].nsel;
			DCFS[eq].fdist=eqkfm0[NFtot].distance;
			DCFS[eq].index_cat=eqkfm0[NFtot].index_cat;
			DCFS[eq].which_pts=eqkfm0[NFtot].selpoints;
			DCFS[eq].t=eqkfm0[NFtot].t;
			M0=0.0;
			for (int f=0; f<Nfaults[eq1]; f++) M0+=pow(10,1.5*(eqkfm0[NFtot+f].mag+6));
			DCFS[eq].m=(2.0/3.0)*log10(M0)-6;
			DCFS[eq].NF=Nfaults[eq1];
			if (DCFS[eq].nsel!=Nsel){
				if (DCFS[eq].nsel>0){
					free_d3tensor(DCFS[eq].S,1,DCFS[eq].nsel,1,3,1,3);
					free_dvector(DCFS[eq].cmb,1,DCFS[eq].nsel);
				}
				DCFS[eq].nsel=Nsel;
				DCFS[eq].S = d3tensor(1,Nsel,1,3,1,3);
				DCFS[eq].cmb=dvector(1,Nsel);
				for (int i=1; i<=Nsel; i++) DCFS[eq].cmb[i]=0.0;
			}
			NFtot+=Nfaults[eq1];
		}

		*DCFS_out=DCFS;
		if(procId == 0) {
			if (flog){
				fprintf(flog,"DCFS structure set up.\n");
				fflush(flog);
			}
		}
    }

    return(0);
}

int setup_afterslip_evol(double Teq, double t0, double t1, double *Cs, double *ts,
						 int Nfun, struct eqkfm **eqk_aft, double *t_afterslip,
						 int Nas,int Nfaults, int afterslip, int *L, double **times2,
						 double **tevol_afterslip, long *seed) {

//if splines are, eqk_aft is substituted with more densily discretized version (L steps instead of Nas).
//Cs, ts, coefficients of temporal evolution functions (See Savage Parkfield paper). Nfun: no. of such funtions.
//t_afterslip indices: [0,1,2,...Nas-1].
// eq_aft has indices: [0...Nas*Nfaults-1].

	int err=0;
	int L0=0, NFL, i;
	int splines;
	double TAU=200000, dtau=(afterslip==0)? 10000 : 700;	//todo allow to set from outside.
	double Kotau;
	double M0,mu;
	double smallstepstime=10;
	double now, prev, norm, curr;
	double Tendaft=t_afterslip[Nas-1];	//Time to which cumulative afterslip snapshot refers.
	struct eqkfm *eq_aft, *eqkfm_aftsplines;

	eq_aft= *eqk_aft;

	splines=(Nas>1)? 1 : 0;
	//todo find L first, then setup vectors.
	*times2=dvector(0,*L-1);
	*tevol_afterslip=dvector(0,*L-1);

	err+=findtimestepsomori(Teq, t0, fmin(smallstepstime,t1), 0, 183, TAU, 0.3*dtau, 0.6, 0.001, (*times2)+1, &Kotau, L);
	if (smallstepstime<t1) err+=findtimestepsomori(Teq, smallstepstime, t1, 0, 183, TAU, dtau, 0.6, 0.001, (*times2)+*L+1, &Kotau, &L0);
	*L+=L0+1;
	(*times2)[0]=Teq-1e-4;

	// Temporal evolution of afterslip.//
	M0=pow(10,1.5*(eq_aft[(Nas-1)*Nfaults].mag+6.0));
	mu=M0/(eq_aft[(Nas-1)*Nfaults].tot_slip*eq_aft[(Nas-1)*Nfaults].L*eq_aft[(Nas-1)*Nfaults].W*1e12);

	if (afterslip!=0){
		if (splines==0){
			norm=0.0;
			for (int i=0; i<Nfun; i++) norm+=Cs[i]*log(1+(Tendaft-Teq)/ts[i]);
			now=0.0;
			curr=0.0;

			for (int t=1; t<=*L; t++){
				prev=curr;
				now= (*times2)[t]-Teq;
				curr=0.0;
				for (int i=0; i<Nfun; i++) curr+=Cs[i]*log(1+now/ts[i]);
				(*tevol_afterslip)[t-1]= (curr-prev)/norm;
			}
		}
		else{
			//todo double check indexing of eqkfm_aftsplines (NB t_afterslip and afterslip_models has been changed, so indices are [0...Nas-1]).
			//double *t_afterslip_dum=dvector(0,L);
			for (int i=0; i<Nas; i++) t_afterslip[i]-=Teq;	//since functions below start from t=0;
			for (int i=0; i<=*L; i++) (*times2)[i]-=Teq;	//since functions below start from t=0;

			NFL=Nfaults*(*L);
			eqkfm_aftsplines=eqkfm_array(0,Nfaults*(*L+1));
			splines_eqkfm(eq_aft-1, Nas, Nfaults, t_afterslip-1, (*times2)-1, *L, eqkfm_aftsplines-1, seed);

			for (int f=0; f<Nfaults; f++) {
				for (int l=*L-1; l>=0; l--) eqkfm_aftsplines[Nfaults*l+f].tot_slip=0.0;
				copy_eqkfm_all(eqkfm_aftsplines[Nfaults*(*L-1)+f], eqkfm_aftsplines+NFL+f);
				eqkfm_aftsplines[NFL+f].slip_str=dvector(1,eq_aft[f].np_di*eq_aft[f].np_st);
				eqkfm_aftsplines[NFL+f].slip_dip=dvector(1,eq_aft[f].np_di*eq_aft[f].np_st);
				for (int p=1; p<=eq_aft[f].np_di*eq_aft[f].np_st; p++) {
					//copy cumulative value at the end of array (will be needed later to calculate random high freq. slip).
					eqkfm_aftsplines[NFL+f].slip_str[p]=eqkfm_aftsplines[Nfaults*(*L-1)+f].slip_str[p];
					eqkfm_aftsplines[NFL+f].slip_dip[p]=eqkfm_aftsplines[Nfaults*(*L-1)+f].slip_dip[p];
					eqkfm_aftsplines[NFL+f].tot_slip+=sqrt(pow(eqkfm_aftsplines[NFL+f].slip_str[p],2)+pow(eqkfm_aftsplines[NFL+f].slip_dip[p],2));
					for (int l=*L; l>1; l--) {
						i=Nfaults*(l-1)+f;
						eqkfm_aftsplines[i].slip_str[p]-=eqkfm_aftsplines[i-Nfaults].slip_str[p];
						eqkfm_aftsplines[i].slip_dip[p]-=eqkfm_aftsplines[i-Nfaults].slip_dip[p];
						eqkfm_aftsplines[i].tot_slip+=sqrt(pow(eqkfm_aftsplines[i].slip_str[p],2)+pow(eqkfm_aftsplines[i].slip_dip[p],2));
					}
				 }
				for (int l=0; l<NFL+Nfaults; l++){
					eqkfm_aftsplines[l].tot_slip*=(1.0/(eqkfm_aftsplines[l].np_di*eqkfm_aftsplines[l].np_st));
					M0=mu*eqkfm_aftsplines[l].tot_slip*eqkfm_aftsplines[l].L*eqkfm_aftsplines[l].W*1e12;
					eqkfm_aftsplines[l].mag=(2.0/3.0)*log10(M0)-6.0;
				}

			}
			*eqk_aft=eqkfm_aftsplines;
			for (int i=0; i<Nas; i++) t_afterslip[i]+=Teq;	//revert shift from before;
			for (int i=0; i<=*L; i++) (*times2)[i]+=Teq;	//revert shift from before;
		}
	}

	else {
		for (int t=1; t<=*L; t++) (*tevol_afterslip)[t-1]=0;
	}

	return(err!=0);
}
