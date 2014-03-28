/*
 * read_eqkfm.c
 *
 *  Created on: Jan 16, 2013
 *      Author: camcat
 */

#include "read_eqkfm.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "../general/mem_mgmt.h"
#include "../general/setup.h"
#include "../general/struct_conversions.h"
#include "../seis/soumod1.h"
#include "../seis/WellsCoppersmith.h"
#include "../util/nrutil.h"
#include "read_eqkfm_fsp.h"

int eqkfm_addslipmodels(struct eqkfm *eqfm1, struct slipmodels_list all_slipmodels, struct eqkfm **eqfm_comb,
						int **which_events, int N1, int *Ncomb, int **nfout,
						double dt, double dmag, double res, struct crust crst,
						int refine, int taper) {
/* assumes that eqfm1 only has single fault events, but slip models may have more...
 * does not do spatial selection (unlike combine_eqkfm).
 * return combined catalog, rather than list of indices mapping one catalog to the other.
 * return error code (1 if not all events from second catalog are selected).
 */
	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int *which_slipmod;
	int N2, N3=0, c_evfound=0;
	int *nf2, nsm, nfaults;
	int c2=0, c3=0;	//counters.
	int err=0, j;
	double d_close=10.0;	//todo read somewhere.
	static struct set_of_models dummy_parentsetofmodels;
	dummy_parentsetofmodels.Nmod=0;	//this value indicates than no slip model is available (will use synthetic slip model from foc. mec. or isotropic field).

	if(procId == 0) {
		if (flog) fprintf(flog, "\nBuilding mainshock slip models (eqkfm_addslipmodels).\n");
	}

	N2=all_slipmodels.NSM;
	nf2=all_slipmodels.Nfaults;

	if (refine==0) taper=0;
	*Ncomb=0;

	//find no. of faults in each slip model, and save the largest value for each event.
	nsm=0;
	for (int i=0; i<all_slipmodels.NSM; i++){
		nfaults=0;
		for (int n=1; n<=all_slipmodels.no_slipmodels[i]; n++){
			if (!(strcmp(cmb_format,"farfalle"))) err+=read_farfalle_eqkfm(all_slipmodels.slipmodels[nsm], NULL, all_slipmodels.Nfaults+i);
			else {
				if (!(strcmp(cmb_format,"pscmp"))) err+=read_pscmp_eqkfm(all_slipmodels.slipmodels[nsm], NULL, all_slipmodels.Nfaults+i);
				else {
					if (!(strcmp(cmb_format,"fsp"))) err+=read_fsp_eqkfm(all_slipmodels.slipmodels[nsm], NULL, all_slipmodels.Nfaults+i);
					else {
						if (flog) fprintf(flog,"Unknown slip model format %s (eqkfm_addslipmodels).\n", cmb_format);
						err=1;
					}
				}
			}
			nsm+=1;
			nfaults=fmax(nfaults,all_slipmodels.Nfaults[i]);
		}
		all_slipmodels.Nfaults[i]=nfaults;
	}

	if (err){
		if(procId == 0) {
			if (verbose_level) printf("Error in reading input files. Exiting.\n");
			if (flog) fprintf(flog, "Error in reading input files. Exiting.\n");
		}

		return 1;
	}

	//todo make sure mmain is already assigned here!! (NB: all this may not be needed if clear association b/w cat and foc mec is given).
	which_slipmod = combine_cats(all_slipmodels.tmain, timesfromeqkfm(eqfm1, N1, (int *) 0),
								 all_slipmodels.mmain, magssfromeqkfm(eqfm1, N1, (int *) 0 ),
								 all_slipmodels.NSM, N1, dt, dmag);

	for (int i=0; i<N1; i++){
		N3= (which_slipmod[i]==-1)? N3+1 : N3+ nf2[which_slipmod[i]];
	}

	*eqfm_comb=eqkfm_array(0,N3-1);
	if (which_events!=NULL) {
		*which_events=ivector(0,N3-1);
	}
	if (nfout) *nfout=ivector(0,N1-1);

	for(int i=0; i<N1; i++) {
		if(eqfm1[i].is_mainshock) {
			eqfm1[i].nsel=0;	//deactivate it as source aftershock;
			(*eqfm_comb)[c3].nsel=crst.N_allP;

			if(which_slipmod[i]==-1) {
				copy_eqkfm_all(eqfm1[i], (*eqfm_comb)+c3);
				if (which_events!=NULL) (*which_events)[c3]=i;
				(*eqfm_comb)[c3].parent_set_of_models=&dummy_parentsetofmodels;
				if (!eqfm1[i].is_slipmodel){
					if(procId == 0) {
						if (verbose_level>0) printf(" ** Warning: slip model or focal mechanism not available for large event at t=%.2lf, mag=%.2lf -> will use isotropic field. **\n", eqfm1[i].t, eqfm1[i].mag);
						if (flog) fprintf(flog, "Warning: slip model or focal mechanism not available for large event at t=%.2lf, mag=%.2lf -> will use isotropic field.\n", eqfm1[i].t, eqfm1[i].mag);
					}
				}
				else {
					err = focmec2slipmodel(crst, (*eqfm_comb)+c3, res, refine, taper);
					if(procId == 0) {
						if (err && verbose_level>0) printf("Error in creating slip model (function: eqkfm_addslipmodels)\n");
						else {
							if (verbose_level>0) printf("Using synthetic slip model from focal mechanism for large event at t=%.5e, mag=%.2lf\n", eqfm1[i].t, eqfm1[i].mag);
							if (flog) fprintf(flog, "Using synthetic slip model from focal mechanism for large event at t=%.5e, mag=%.2lf\n", eqfm1[i].t, eqfm1[i].mag);
						}
					}
				}
				(*nfout)[*Ncomb]=1;
				c3+=1;
				*Ncomb+=1;
			}
			else {
				j=which_slipmod[i];
				nsm=0;
				for (int n=0; n<j; n++) nsm+=all_slipmodels.no_slipmodels[n];
				c2=0;
				if(procId == 0) {
					if (verbose_level>0) printf("Using slip model %s from focal mechanism for large event at t=%.5e, mag=%.2lf\n", all_slipmodels.slipmodels[nsm], eqfm1[i].t, eqfm1[i].mag);
					if (flog) fprintf(flog, "Using slip model %s from focal mechanism for large event at t=%.5e, mag=%.2lf\n", all_slipmodels.slipmodels[nsm], eqfm1[i].t, eqfm1[i].mag);
				}

				err += setup_eqkfm_element((*eqfm_comb)+c3, all_slipmodels.slipmodels+nsm,
										   all_slipmodels.no_slipmodels[j], crst.mu, all_slipmodels.disc[j],
										   all_slipmodels.tmain[j], d_close, crst.N_allP, crst.list_allP,
										   all_slipmodels.mmain+j, 1, 1, NULL, crst.lat0, crst.lon0);

				if (which_events!=NULL) (*which_events)[c3]=i;
				if (nfout)(*nfout)[*Ncomb]=nf2[which_slipmod[i]];
				for (int cc3=c3; cc3<c3+nf2[which_slipmod[i]]; cc3++) (*eqfm_comb)[c3].distance=eqfm1[i].distance;
				c3+=nf2[which_slipmod[i]];
				c_evfound+=1;
				*Ncomb+=1;
			}
		}
	}

	if (c_evfound<N2){
		if(procId == 0) {
			if (verbose_level>0) printf("** Warning: some slip models are not used, since events were not found in catalog (or are outside time boundaries)! (function: eqkfm_addslipmodels)\n");
			if (flog) fprintf(flog, "Warning: some slip models are not used, since events were not found in catalog (or are outside time boundaries)! (function: eqkfm_addslipmodels)\n");
		}
	}

	return (err);
}

int focmec2slipmodel(struct crust crst, struct eqkfm *eqfm1, double res, int refine, int taper){
	//eqfm1 should already contain strX, dipX, rakeX, mag, whichfm.

	struct eqkfm eqkfm0;
	struct eqkfm *eqkfmP;
	int err=0;
	double slip;

	if (refine) {
		copy_eqkfm_all((*eqfm1), &eqkfm0);
		eqkfmP=&eqkfm0;
	}
	else eqkfmP=eqfm1;

	(*eqkfmP).is_slipmodel=1;
	(*eqkfmP).np_st=1;
	(*eqkfmP).np_di=1;
	(*eqkfmP).slip_str=dvector(1,1);
	(*eqkfmP).slip_dip=dvector(1,1);
	//(*eqkfmP).taper=dvector(1,4);
	//	for (int i=1; i<=4; i++) (*eqkfmP).taper=1.0; //probably not used (aftershocks are not resampled or tapered).
	(*eqkfmP).pos_s=dvector(1,1);	//location of patches within fault; [0], [0] for single patch events.
	(*eqkfmP).pos_d=dvector(1,1);
	(*eqkfmP).pos_s[1]=0;	//location of patches within fault; [0], [0] for single patch events.
	(*eqkfmP).pos_d[1]=0;

	//TODO in theory, should find L and W for both mech. (rake differs).
	WellsCoppersmith((*eqkfmP).mag, (*eqkfmP).rake1, &((*eqkfmP).L), &((*eqkfmP).W), &slip);
	slip=(*eqkfmP).tot_slip=pow(10,(1.5*((*eqkfmP).mag+6)))*(1.0/(crst.mu*pow(10,12)*(*eqkfmP).W*(*eqkfmP).L));
	switch ((*eqkfmP).whichfm){
		case 1:
			if ((*eqkfmP).depth<0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1)) {
				if (verbose_level>1) printf("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
				(*eqkfmP).depth=0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1);
			}
			(*eqkfmP).slip_str[1]=slip*cos(DEG2RAD*(*eqkfmP).rake1);
			(*eqkfmP).slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfmP).rake1);
			break;
		case 2:
			if ((*eqkfmP).depth<0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip2)) {
				if (verbose_level>1) printf("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
				(*eqkfmP).depth=0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip2);
			}
			(*eqkfmP).slip_str[1]=slip*cos(DEG2RAD*(*eqkfmP).rake2);
			(*eqkfmP).slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfmP).rake2);
			break;
		case 0:
			if ((*eqkfmP).depth<0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1)) {
				if (verbose_level>1) printf("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
				(*eqkfmP).depth=0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1);
			}
			if (verbose_level>1) printf("Warning: ambiguous nodal plane (focmec2slipmodel).\n");
			(*eqkfmP).slip_str[1]=slip*cos(DEG2RAD*(*eqkfmP).rake1);
			(*eqkfmP).slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfmP).rake1);
			break;
		default:
			if (verbose_level>0) printf("Error: illegal value for whichfm (focmec2slipmodel). Exiting. \n");
			return(1);
	}

	if (refine) {
		err+=suomod1_resample(eqkfm0, eqfm1, res, 0.0);	//create a slip model with right resolution.
		if (taper) err+=suomod1_taper((*eqfm1), eqfm1);
	}

	return err;
}


//--------------------------read slip model files----------------------//

int read_eqkfm(char *fname, struct eqkfm **eqfm1, int *NF_out, double *Mw, double mu) {
	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int NF, NP, err=0;
	double slip, dip_slip=0.0, strike_slip=0.0, M0=0.0;

	if (!(strcmp(cmb_format,"farfalle"))) err=read_farfalle_eqkfm(fname, eqfm1, &NF);
	else {
		if (!(strcmp(cmb_format,"pscmp"))) err=read_pscmp_eqkfm(fname, eqfm1, &NF);
		else{
			if (!(strcmp(cmb_format,"fsp"))) err+=read_fsp_eqkfm(fname, eqfm1, &NF);
			else {
				if(procId == 0) {
					if (flog) fprintf(flog,"Unknown slip model format %s (read_eqkfm).\n", cmb_format);
					if (verbose_level>0) printf(" ** Error: could not read slip model file %s (read_slipmodel). ** \n ", fname);
				}
				err=1;
			}
		}
	}

	if (err==1){
		if(procId == 0) {
			if (verbose_level>0) printf(" ** Error: could not read slip model file %s (read_slipmodel). ** \n ", fname);
		}

		return (1);
	}

	if (eqfm1){
		for (int f=0; f<NF; f++){
			(*eqfm1)[f].str2=(*eqfm1)[f].str1;
			(*eqfm1)[f].dip2=(*eqfm1)[f].dip1;	//in this case there is no ambiguity (correct plane is known).
			(*eqfm1)[f].whichfm=1;
			(*eqfm1)[f].is_slipmodel=1;
			(*eqfm1)[f].is_mainshock=1;
			(*eqfm1)[f].noise=0;
			(*eqfm1)[f].taper=ivector(1,4);
			(*eqfm1)[f].tot_slip=0.0;
			NP=(*eqfm1)[f].np_di*(*eqfm1)[f].np_st;
			for (int p=1; p<=NP; p++) {
				strike_slip+=(*eqfm1)[f].slip_str[p];
				dip_slip+=(*eqfm1)[f].slip_dip[p];
				slip=pow((*eqfm1)[f].slip_str[p]*(*eqfm1)[f].slip_str[p]+(*eqfm1)[f].slip_dip[p]*(*eqfm1)[f].slip_dip[p],0.5);
				(*eqfm1)[f].tot_slip+=slip/NP;
			}
			(*eqfm1)[f].rake1=(*eqfm1)[f].rake2=(-RAD2DEG)*atan2(dip_slip, strike_slip);
			if (Mw) M0+=mu*(*eqfm1)[f].tot_slip*(*eqfm1)[f].L*(*eqfm1)[f].W*1e12;
			(*eqfm1)[f].mag=(2.0/3.0)*log10(mu*(*eqfm1)[f].tot_slip*(*eqfm1)[f].L*(*eqfm1)[f].W*1e12)-6.0;
		}
	}
	if (Mw) {
		if (!eqfm1){
			if(procId == 0) {
				if (verbose_level) printf("*Warning: can not calculate Mw since eqfm1==NULL (read_eqkfm).*\n");
				if (flog) printf("Warning: can not calculate Mw since eqfm1==NULL (read_eqkfm).\n");
			}
		}
		else *Mw=(2.0/3.0)*log10(M0)-6.0;
	}
	if (NF_out) *NF_out=NF;
	return (0);
}

int read_farfalle_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF_out) {
	// [Fahad] Variables used for MPI.
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	FILE *fin;
	struct eqkfm *eqfm;
	int NF, NP, nchar=200;
	char line[nchar];
	double dlen, dwid;
	double *slips, *rakes;
	double plength, pwidth;
	long file_pos;

	if(procId == 0) {
		fin = fopen(fname,"r");
		if(fin == NULL) {
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		if(procId == 0) {
			if (verbose_level>0) printf("Invalid input file (%s).\n", fname);
			if (flog) fprintf(flog,"Invalid input file passed to read_farfalle_eqkfm (%s).\n", fname);
		}

		return (1);
	}
	else {
		if(procId == 0) {
			fgets(line,nchar,fin);
			fgets(line,nchar,fin);
			sscanf(line,"%d",&NF);
		}

		#ifdef _CRS_MPI
			MPI_Bcast(&NF, 1, MPI_INT, 0, MPI_COMM_WORLD);
		#endif

			if (eqfm_out) {
			eqfm=eqkfm_array(0,NF-1);
			for (int f=0; f<NF; f++){
				if(procId == 0) {
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%lf %lf %lf",&(eqfm[f].lat), &(eqfm[f].lon), &(eqfm[f].depth));
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
	//				eqfm[f].whichfm=1;	// [Fahad] Move it after the following statements for MPI
					sscanf(line,"%lf %lf",&(eqfm[f].str1), &(eqfm[f].dip1));
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%lf %lf",&(eqfm[f].L), &plength);
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%lf %lf",&(eqfm[f].W), &pwidth);
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%d %d",&(eqfm[f].np_st), &(eqfm[f].np_di));
				}
				eqfm[f].whichfm=1;
			#ifdef _CRS_MPI
				MPI_Bcast(&(eqfm[f].lat), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].lon), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].depth), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].str1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].dip1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].L), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&plength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].W), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&pwidth, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].np_st), 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].np_di), 1, MPI_INT, 0, MPI_COMM_WORLD);
			#endif
				NP=eqfm[f].np_st*eqfm[f].np_di;
				eqfm[f].pos_s=dvector(1,NP);
				eqfm[f].pos_d=dvector(1,NP);
				eqfm[f].slip_str=dvector(1,NP);
				eqfm[f].slip_dip=dvector(1,NP);
				if ((f==0) | (NP>eqfm[f-1].np_st*eqfm[f-1].np_di)){
					if (f!=0){
						free_dvector(slips, 1,1);
						free_dvector(rakes, 1,1);
					}
					slips=dvector(1,NP);
					rakes=dvector(1,NP);
				}
				if(procId == 0) {
					fgets(line,nchar,fin);
					file_pos=ftell(fin);
				}
				dlen=eqfm[f].L/eqfm[f].np_st;
				dwid=eqfm[f].W/eqfm[f].np_di;
				for (int w=1; w<=eqfm[f].np_di; w++){
					for (int l=1; l<=eqfm[f].np_st; l++){
						if(procId == 0) {
							fscanf(fin, "%lf", slips+(w-1)*eqfm[f].np_st+l);
						}
											//if points refer to center of the faut:
					//eqfm[f].pos_s[(w-1)*eqfm[f].np_st+l]=-0.5*eqfm[f].L+(l-0.5)*dlen;
					//eqfm[f].pos_d[(w-1)*eqfm[f].np_st+l]=-0.5*eqfm[f].W+(w-0.5)*dwid;

					//if points refer to top left corner of the fault:
					//eqfm[f].pos_s[(w-1)*eqfm[f].np_st+l]=(l-0.5)*dlen;
					//eqfm[f].pos_d[(w-1)*eqfm[f].np_st+l]=(w-0.5)*dwid;

					/* Abi's email "The hypocenter is the reference point. Then, you draw a line in the strike direction,
					 * and the distance between the hypocenter and the end of the fault is the partial length.
					 * The partial width is the distance to the top of the fault.
					 * It is the end that goes in the opposite direction of the strike, starting from the reference point, if I remember well,
					 * so that you always start drawing the fault in the direction of the strike.*/
					eqfm[f].pos_s[(w-1)*eqfm[f].np_st+l]=(l-0.5)*dlen-plength;
					eqfm[f].pos_d[(w-1)*eqfm[f].np_st+l]=(w-0.5)*dwid-pwidth;
					}
				}

				#ifdef _CRS_MPI
					MPI_Bcast(slips, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				#endif

				if(procId == 0) {
					fseek(fin,file_pos, SEEK_SET);		//since fscanf sometimes reads until next line, sometimes not (depending on small differences in the file, e.g. spaces).
					for (int w=1; w<=eqfm[f].np_di; w++) fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					for (int w=1; w<=eqfm[f].np_di; w++) fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					file_pos=ftell(fin);
					for (int w=1; w<=eqfm[f].np_di; w++){
						for (int l=1; l<=eqfm[f].np_st; l++){
							fscanf(fin, "%lf", rakes+(w-1)*eqfm[f].np_st+l);
						}
					}
					fseek(fin,file_pos, SEEK_SET);		//since fscanf sometimes reads until next line, sometimes not (depending on small differences in the file, e.g. spaces).
					for (int w=1; w<=eqfm[f].np_di; w++) fgets(line,nchar,fin);
				}

				#ifdef _CRS_MPI
					MPI_Bcast(rakes, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				#endif

				for (int p=1; p<=NP; p++) {
					eqfm[f].slip_str[p]=slips[p]*cos(DEG2RAD*rakes[p]);
					eqfm[f].slip_dip[p]=-slips[p]*sin(DEG2RAD*rakes[p]);
				}
			}

			*eqfm_out=eqfm;
		}
		*NF_out=NF;
	}

	return(0);
}

int read_pscmp_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF2){
//mu in same units as in crst.

	// [Fahad] Variables used for MPI.
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	FILE *fin;
	int dumerror;
	double junk;
	int NP, NF, djunk;
	int nchar=500;
	char cjunk[100];
	struct eqkfm *eqfm1;
	char comm[]="#";
	char line[nchar];

	if(procId == 0) {
		fin = fopen(fname,"r");
		if(fin == NULL) {
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		if(procId == 0) {
			if (verbose_level>0) printf("Invalid input file (%s).\n", fname);
			if (flog) fprintf(flog,"Invalid input file passed to read_pscmp_eqkfm (%s).\n", fname);
		}

		return (1);
	}

	line[0]=comm[0];
	if(procId == 0) {
		while (line[0]==comm[0])fgets(line, nchar, fin);
		for (int i=1; i<=3; i++){
			while (line[0]!=comm[0])fgets(line, nchar, fin);
			while (line[0]==comm[0])fgets(line, nchar, fin);
		}
		sscanf(line, "%d", &NF);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NF, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (eqfm_out) {
		eqfm1=eqkfm_array(0,NF-1);
		line[0]=comm[0];
		if(procId == 0) {
			while (line[0]==comm[0])fgets(line, nchar, fin);
		}
		for (int f=0; f<NF; f++){
			if(procId == 0) {
				sscanf(line, "%d   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %d   %d   %lf",
						&djunk, &(eqfm1[f].lat), &(eqfm1[f].lon), &(eqfm1[f].depth), &(eqfm1[f].L),
						&(eqfm1[f].W), &(eqfm1[f].str1), &(eqfm1[f].dip1), &(eqfm1[f].np_st),
						&(eqfm1[f].np_di), &(eqfm1[f].t));
				fgets(line, nchar, fin);
			}

			#ifdef _CRS_MPI
				MPI_Bcast(&(eqfm1[f].lat), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].lon), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].depth), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].L), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].W), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].str1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].dip1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].np_st), 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].np_di), 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].t), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			#endif

			eqfm1[f].str2=eqfm1[f].str1; eqfm1[f].dip2=eqfm1[f].dip1;	//in this case there is no ambiguity (correct plane is known).
			eqfm1[f].whichfm=1;
			NP=eqfm1[f].np_st*eqfm1[f].np_di;
			eqfm1[f].pos_s=dvector(1,NP);
			eqfm1[f].pos_d=dvector(1,NP);
			eqfm1[f].slip_str=dvector(1,NP);
			eqfm1[f].slip_dip=dvector(1,NP);
			if(procId == 0) {
				for (int p=1; p<=NP; p++) {
						sscanf(line, "%lf   %lf    %lf   %lf   %lf",
								&(eqfm1[f].pos_s[p]), &(eqfm1[f].pos_d[p]),
								&(eqfm1[f].slip_str[p]), &(eqfm1[f].slip_dip[p]),
								&junk);
						if (!feof(fin)) fgets(line, nchar, fin);
				}
			}

			#ifdef _CRS_MPI
				MPI_Bcast(eqfm1[f].pos_s, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(eqfm1[f].pos_d, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(eqfm1[f].slip_str, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(eqfm1[f].slip_dip, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			#endif
		}
	}

	if(procId == 0) {
		fclose(fin);
	}

	if (NF2!= NULL) *NF2=NF;
	if (eqfm_out) *eqfm_out=eqfm1;

	return (0);
}
