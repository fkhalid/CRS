/*
 * main.c
 *
 *  Created on: Aug 23, 2013
 *      Author: camcat
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include "src/defines.h"
#include "src/seis/background_rate.h"
#include "src/general/CRS_LogLikelihood.h"
#include "src/inp_out/print_output.h"
#include "src/inp_out/propagate_results.h"
#include "src/inp_out/read_crust.h"
#include "src/inp_out/read_csep_template.h"
#include "src/inp_out/read_eqkfm.h"
#include "src/inp_out/read_inputfile.h"
#include "src/seis/GR.h"
#include "src/util/hash.h"
#include "src/util/moreutil.h"
#include "src/util/nrutil.h"

//todo make sure first event is selected (criterion: t>=tstart, not t>tstart!!).
//todo make sure output comments make sense at various verbosity levels.
//todo introduce tapering also when there are no uncertainties? also make sure mainshocks are tapered.
//todo make program work for non uniform grid (set griderr=0 in that case).
//todo check that lat=(-180, 180) are interpreted correctly throughout the program.

int verbose_level=1;	//todo allow passing as argument.
double DCFS_cap;
int gridPMax;
char cmb_format[120];
int CSEPmode=0;
FILE *flog=NULL;

void error_quit(char * message);

int main (int argc, char **argv){

	static int clean=0;

	setenv("TZ", "UTC", 1);
	int run_tests=0;

	if (run_tests){
		verbose_level=2;
		test_allOkada_simple_multiplerec();
		printf("Done!\n");
		return (0);
	}

	int err=0;

	FILE *fout, *fin;

	int LLinversion, forecast, extra_output;
	int Nchar=120, Nchar_long=500;
	char focmec_format[Nchar];
	char fname[Nchar], msg[Nchar],	infile[Nchar], crust_file[Nchar], fore_template[Nchar], catname[Nchar], outname[Nchar],
		syscopy[Nchar], background_rate_file[Nchar], slipmodelfile[Nchar], afterslipmodelfile[Nchar], modelparametersfile[Nchar],
		logfile[Nchar], print_LL[Nchar], outnamemod[Nchar], error_msg[120];
	char print_cmb[Nchar],  print_forex[Nchar],  print_foret[Nchar],  printall_cmb[Nchar],  printall_forex[Nchar],  printall_foret[Nchar];
	char line[Nchar_long];
	char **focmeccats;

	//Slip model + catalog variables:
    struct slipmodels_list all_slipmodels, all_aslipmodels;
	struct eqkfm *eqkfm1=0, *eqkfm0res=0, *eqkfm_aft=0;	//contain respectively: aftershocks with fm, aftershocks, mainshocks, mainshocks resampled, afterslip, afterslip resampled.
	double Mag_main;
	double dt, dM, xytoll, ztoll, border;
	double Hurst;

	int Nm, Nas=0, *Nfaults_all=0;
	double res, gridresxy, gridresz;
	struct catalog cat, cat2;
	int Ntot;
	double **focmec=0;
	int *fmzonelimits;

	//DCFS variables:
	struct pscmp *DCFS;
	double dDCFS;
	struct crust crst;
	int NgridT, NFM=0;	//number of gridpoints, number of aftershocks for which foc.mec is known, number of foc.mec from past seismicity.
	struct Coeff_LinkList *AllCoeff;
    int *which_main;
    int sample_allFM=0, aftershock_mode, aftershocksMC, load_focmec;
    double  Mc_source;

	//RateState variables
	int fixta, fixAsig, fixr;
	double Asig_min, Asig_max, Asig0, min_dAsig;
	double ta_min, ta_max, ta0, min_dta;
	double r0;
	int nAsig0, nAsig;
	int nta0, nta;
	double *maxAsig, *maxta, *maxr;
	double Asig, dAsig, ta, dta, r;

	struct flags flags;
	int Nsur, Nslipmod;
	int no_fm_cats;
	//int sample_allFM=0;

	struct tm reftime, times;
	double tstartLL, tendLL, Tend, Tstart, tw, t_oldsnap, tstart_calc;
	double fore_dt;
	double *tts;
	int Ntts;
	double extra_time;
	double t_firstmain;
	double smoothing;	//for calculating background seismicity.
	double t_back;		//for calculating background seismicity.

	//grid search variables:
	int p=0, p_found;
	int snapshot_exists, use_snap, multi_gammas, use_bg_rate, use_bg_rate_file;
	long old_hash;
	double *LLs, **LLs_old_matrix, *Ldums0, *Nev, *I, LLmax;
	double *gamma_bgrate=NULL;
	double 	**gammas=NULL, \
			**gammas_new=NULL, \
			**gammas_old=NULL, \
			**gammas_maxLL=NULL;

	//hash (to propagate results in CSEP mode).
	char *str_data;
	int npar_bool, npar_int, npar_double, str_data_l;
	long new_hash;

	//to switch between slip models.
	int refresh;
	int slipmodel_combinations=1;
	int *dim;
	int *Nsm, nf=0;

	//temporary variables.
	double minmag;
	long seed;
	int nc;
	int first_main;
	double tnow;
	int j0, N, N_min_events, current_ev;

	printf("Input file: %s\n", argv[1]);
	sscanf(argv[1],"%s", infile);

	if (strcmp(infile,"clean")==0) clean=1;

	if (CSEPmode && clean){
		snapshot_exists=check_if_snapshot_exists(".", NULL, reftime, NULL);
		if (snapshot_exists) {
			remove(check_if_snapshot_filename);
			printf("removed file %s\n",check_if_snapshot_filename);
			if (CSEPmode) {
				sprintf(fname,"rm %s/*",old_LLfolder);	//platform specific...
				system(fname);
				printf("%s\n",fname);
			}
		}
		return 0;
	}

	//-----------------------read input file -------------------//

	err=read_inputfile(infile, outname, NULL, crust_file, fore_template, catname, &focmeccats, background_rate_file,
			slipmodelfile, afterslipmodelfile,	modelparametersfile, logfile, &extra_output, &reftime, &Tstart, &Tend, &seed,
			cmb_format, &no_fm_cats);

	if (err) {
		sprintf(msg,"Error reading input file %s.\n", infile);
		error_quit(msg);
	}
	//set things for Csep mode:
	if (CSEPmode) {
		sprintf(focmec_format, "CSEP");
		sprintf(cmb_format, "farfalle");
	}
	else {	//todo only leave CSEP format for final version.
		nc=countcol(focmeccats[0]);
		if ((nc==7) | (nc==10)) sprintf(focmec_format,"7col");
		else sprintf(focmec_format,"CSEP");
	}

//-----------------------read model parameters-------------------//

	read_modelparmeters(modelparametersfile, reftime, &N_min_events, &fixr, &fixAsig, &fixta, &r0, &Asig0, &ta0,
			&Asig_min, &Asig_max, &ta_min, &ta_max, &nAsig0, &nta0, &tstartLL, &extra_time, &tw, &fore_dt, &t_back,
			&Nsur, &Nslipmod, &flags, &Hurst, &Mc_source, &use_bg_rate, &(cat.Mc), &Mag_main, &DCFS_cap, &gridPMax,
			&dt, &dM, &xytoll, &ztoll, &border, &res, &gridresxy, &gridresz, &smoothing, &LLinversion, &forecast);

//------- change flags if input files are missing:----//

	if (!focmeccats && flags.err_recfault){
		if (verbose_level>0) printf("Warning: InputCatalogFocMecFile not given: will not use variable receiver faults.\n");
		flags.err_recfault=0;
	}
	if ((strcmp(afterslipmodelfile,"")==0) && flags.afterslip){
		if (verbose_level>0) printf("Warning: InputListAfterslipModels not given: will not use afterslip.\n");
		flags.afterslip=0;
	}
	if (strcmp(background_rate_file,"")==0)	use_bg_rate_file=0;
	else use_bg_rate_file=1;

//----------------------------------------------------------------------------------------------//

	if (strcmp(logfile,"")!=0){
		sprintf(syscopy,"date > %s", logfile);
		system(syscopy);
		flog=fopen(logfile,"a");
		fprintf(flog,"!--------------------------!\n!------input file: --------!\n!--------------------------!\n");
		fclose(flog);
		sprintf(syscopy,"cat %s >> %s",infile, logfile);
		system(syscopy);
		flog=fopen(logfile,"a");
		fprintf(flog,"\n\n!--------------------------!\n!------param file: --------!\n!--------------------------!\n");
		fclose(flog);
		sprintf(syscopy,"cat %s >> %s", modelparametersfile, logfile);
		system(syscopy);
		flog=fopen(logfile,"a");
		fprintf(flog,"\nRunning in CSEP mode...\n");
		fflush(flog);
	}
	else flog=NULL;

	sprintf(syscopy,"cp %s %s_inputfile.txt",infile, outname);
	system(syscopy);
	sprintf(syscopy,"cp %s %s_parameters.txt", modelparametersfile, outname);
	system(syscopy);

	err=read_crust(crust_file, fore_template, &crst, gridresxy, gridresz);
	if (err) {
		sprintf(error_msg, "Errors while reading crust file %s or template file %s. Exiting.", crust_file, fore_template);
		error_quit(error_msg);
	}
	if (flags.err_recfault) read_fmindex(crst, fore_template, &(crst.fmzone), &(crst.nofmzones));
	else crst.nofmzones=1;
	NgridT=crst.N_allP;

	dDCFS=(fixAsig)? 0.01*Asig0 : 0.01*Asig_min;	//minimum stress for which points are considered. todo test this (may need smaller value if aftershocks have important cumulative effect).
	if (fixAsig) nAsig0=0;
	if (fixta) nta0=0;

	if (flog) fprintf(flog, "dDCFS (min value for which calculation is done) = %.2e Pa\n", dDCFS);


//---------------------------------------------//
//--------------Setup afterslip----------------//
//---------------------------------------------//

	if (flags.afterslip !=0) {
		read_listslipmodel(afterslipmodelfile, reftime, &all_aslipmodels, res, 1);
		err=setup_eqkfm(all_aslipmodels, crst, flags.err_slipmodel, &eqkfm_aft);
		if (err!=0) error_quit("Error in setting up afterslip slip model - exiting.");
	}
	else eqkfm_aft=NULL;

	flags.splines= (flags.afterslip)? (Nas>1): 0;

//----------------------------------------------------------//
//--------------Setup aftershocks, mainshocks --------------//
//----------------------------------------------------------//

	read_listslipmodel(slipmodelfile, reftime, &all_slipmodels, res, 0);

	aftershocksMC = (flags.aftershocks && flags.full_field && !flags.aftershocks_fixedmec);
	load_focmec= (flags.err_recfault || aftershocksMC);

	if (load_focmec)
		 err=setup_catalogetc(catname, focmeccats, no_fm_cats, focmec_format, reftime, dDCFS, Mag_main, crst, &cat, &eqkfm1, &focmec, &fmzonelimits, flags, &NFM, &Ntot, &Nm, dt, dM,  xytoll, ztoll, border, tw, tstartLL, tstartLL-extra_time, 30);
	else err=setup_catalogetc(catname, focmeccats, no_fm_cats, focmec_format, reftime, dDCFS, Mag_main, crst, &cat, &eqkfm1,   NULL , NULL, flags, NULL, &Ntot, &Nm, dt, dM,  xytoll, ztoll, border, tw, tstartLL, tstartLL-extra_time, 30);
	if (err!=0) error_quit("**Error in setting up catalog or associating events with mainshocks. Exiting. **");

	if (flags.err_recfault && (no_fm_cats!=crst.nofmzones)){
		if (crst.nofmzones>no_fm_cats){
			if (verbose_level) printf("**Error: not enough catalogs of focal mechanisms given! (%d given, at least %d required)**\n", no_fm_cats, crst.nofmzones);
			if (flog) fprintf(flog, "**Error: not enough catalogs of focal mechanisms given! (%d given, at least %d required)**\n", no_fm_cats, crst.nofmzones);
			crst.nofmzones=1;
		}
		else {
			if (verbose_level>1) printf("**Warning: some catalogs of focal mechanisms not used! (%d given, %d used)**\n", no_fm_cats, crst.nofmzones);
			if (flog) fprintf(flog, "**Warning: some catalogs of focal mechanisms not used! (%d given, %d used)**\n", no_fm_cats, crst.nofmzones);
		}
	}

//----------------------------------------------------------//
//-----------------Setup LL inversion period ---------------//
//----------------------------------------------------------//

	first_main=j0=1;
	tnow=-1.0;
	current_ev=N=0;

	while (current_ev<Ntot && tnow<0 && N<N_min_events){
		if (eqkfm1[current_ev].is_mainshock) {
			if (first_main)	{
				tnow=t_firstmain=eqkfm1[current_ev].t;
				first_main=0;
			}
			if (tnow<eqkfm1[current_ev].t){
				for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<eqkfm1[current_ev].t) N+=1;
				j0+=N;
			}
			tnow=eqkfm1[current_ev].t+tw;
		}
		current_ev+=1;
	}
	if (tnow<0){
		for(int j=j0;j<=cat.Z;j++) if(cat.t[j]>=tnow && cat.t[j]<0) N+=1;
	}

	if (N>N_min_events){
		if (flog) fprintf(flog, "\n%d events from catalog can be used for LL inversion - enough to perform inversion (Nmin=%d).\n", N, N_min_events);
		tstartLL=t_firstmain;
		tendLL=0;
	}

	else{
		if (flog) fprintf(flog, "\n%d events from catalog can be used for LL inversion - not enough to perform inversion (Nmin=%d).\n Searching past mainshocks...\n", N, N_min_events);
		if (verbose_level>0) printf("** Warning: fewer than %d events in catalog can be used for LL inversion: will try to find large in the past to fit parameters... **\n", N_min_events);
		cat2.Mc=Mag_main;
		err+=readZMAP(&cat2, NULL, NULL, catname, crst, reftime, t_back, tstartLL, t_back, tstartLL, Mag_main, 0, 0, 0, dDCFS, 1);
		if (cat2.Z!=0){
			if (flog) fprintf(flog, "\nMainshocks found. Reloading catalog for new time period...\n");
			tstartLL=cat2.t[cat2.Z];
			tendLL=fmin(0.0, tstartLL+200);
			free_cat(cat);
			free_eqkfmarray(eqkfm1, 0, Ntot-1);
			err=setup_catalogetc(catname, focmeccats, no_fm_cats, focmec_format, reftime, dDCFS, Mag_main, crst, &cat, &eqkfm1, NULL , NULL, flags, NULL, &Ntot, &Nm, dt, dM,  xytoll, ztoll, border, tw, tstartLL, tstartLL-extra_time, 30);
			if (err!=0) error_quit("**Error in setting up catalog or associating events with mainshocks. Exiting. **");
			free_cat(cat2);
		}
		else {
			if (verbose_level>0) printf("** Warning: No mainshocks found in catalog: will not do parameter estimation, but use default parameters.**\n");
			if (flog) {
				fprintf(flog, "\nMainshocks not found found. Will not perform LL parameter inversion, but use default values: Asig=%.2lf, ta=%.2lf\n", Asig0, ta0);
				fflush(flog);
			}
			LLinversion=0;
		}
	}

	if (load_focmec){
		select_fm_time(focmec, &NFM, Tstart);
		if (!NFM) {
			if (flog) {
				fprintf(flog,"\nNo focal mechanisms available before t=%.2lf (ForecastStartDate). Will not use multiple receiver faults.\n", Tstart);
				if (flags.aftershocks && flags.full_field && !flags.aftershocks_fixedmec)
					fprintf(flog,"No focal mechanisms available before t=%.2lf (ForecastStartDate). Will not use MC sampling of focal planes for aftershochs.\n", Tstart);
				fflush(flog);
			}
			flags.err_recfault=0;
			flags.aftershocks_fixedmec=1;
		}
		else if (flog) {
			fprintf(flog,"\nWill use %d receiver focal mechanisms up to t=%.2lf (ForecastStartDate)\n", NFM,  Tstart);
			fflush(flog);
		}
	}

//----------------------------------------------------------//
//-----------------Add mainshock slip models ---------------//
//----------------------------------------------------------//

	err=eqkfm_addslipmodels(eqkfm1, all_slipmodels, &eqkfm0res, &which_main, Ntot, &Nm, &Nfaults_all, dt, dM, res, crst, 1, 1);
	if (err!=0) error_quit("**Error in setting up catalog or associating events with mainshocks. Exiting. **");

	if (flog && LLinversion){
		fprintf(flog, "Inversion time period: [%2.lf - %2.lf]days, ", tstartLL, tendLL);
		fprintf(flog, "starting with mainshock: t=%.2lf, Mw=%.2lf, lat=%.2lf, lon=%.2lf\n", eqkfm0res[0].t, eqkfm0res[0].mag, eqkfm0res[0].lat, eqkfm0res[0].lon);
		fflush(flog);
	}

	if (flags.aftershocks==0) {
		Ntot=0;
		for (int i=0; i<Nm; i++) which_main[i]=i;	//since DCFS contains only mainshocks.
	}

	//--------------Setup Coefficients and DCFS struct--------------//

	//todo move below (after checking if snapshot exists, and if so avoid Okada calculations).
	setup_CoeffsDCFS(&AllCoeff, &DCFS, crst, eqkfm0res, eqkfm1, Nm, Ntot, Nfaults_all, which_main);

	//--------------------------------------------------------------//
	//					Setup other things							//
	//--------------------------------------------------------------//

	if (!aftershocksMC && !flags.err_recfault && !flags.err_slipmodel && !flags.err_gridpoints) Nsur=Nslipmod=1;	//since there are not sources of uncertainties.
	if (flags.err_recfault && no_fm_cats==1 && !flags.err_slipmodel && Nsur/Nslipmod>NFM) {
	    	Nsur=NFM;
	    	sample_allFM=1;	//To avoid doing more iterations than needed.
	    }

	if (!flags.err_recfault && flags.OOPs) flags.err_recfault=2;	//by convention, this value means OOPs when passed to calculateDCFSrandomized.
	if (!flags.err_recfault) {
		crst.nofmzones=1;
		crst.fmzone=NULL;
	}

	if (!crst.uniform && flags.err_gridpoints){
		if (verbose_level) printf("** Warning: grid is not uniform -> grid error will not be implemented. **\n");
		if (flog) fprintf(flog, "** Warning: grid is not uniform -> grid error will not be implemented. **\n");
		flags.err_gridpoints=0;
	}

	//--------------------------------------------------------------//
	// 					Setup time steps;							//
	//--------------------------------------------------------------//

	if (verbose_level>0) printf("Setting up time steps...");

	//todo should allow this to be set from outside
	double *Cs, *ts;
	int Nfun, L=(flags.afterslip)? 10000 : 2;
	double *times2, *tevol_afterslip=NULL;

	Nfun=1;
	Cs=dvector(0,Nfun-1);
	ts=dvector(0,Nfun-1);
	Cs[0]=436.63;
	ts[0]=14.2653;

	if (flags.afterslip){
		err=setup_afterslip_evol(all_slipmodels.tmain[0], all_slipmodels.tmain[0]+1e-4, fmax(tendLL, Tend), Cs, ts, Nfun, &eqkfm_aft, all_aslipmodels.tmain,
				all_aslipmodels.NSM, all_aslipmodels.Nfaults[0], flags.afterslip, &L, &times2, &tevol_afterslip, &seed);
		if(err) return 1;
	}

	else {
		times2=dvector(0,L);
		times2[0]=fmin(tstartLL, 0.0)-1e-4;
		for (int i=1; i<=L; i++) times2[i]=times2[i-1]+(fmax(tendLL, Tend)+1e-4-times2[0])/(L-1);
	}

	if (verbose_level>0) printf("done\n");
	if (flog) {
		fprintf(flog, "\nSetting up time steps for calculations: %d time steps between times [%.2lf, %.2lf].\n", L, times2[0], times2[L]);
		fflush(flog);
	}

	//--------------------------------------------------------------//
	//			calculate hash corresponding to input files:		//
	//--------------------------------------------------------------//

	npar_bool=10;
	npar_int=2;
	npar_double=1;

	if (CSEPmode){
		str_data_l=npar_bool+3*npar_int+11*npar_double+strlen(catname)+strlen(crust_file)+strlen(fore_template)+strlen(background_rate_file);
		for (int nn=0; nn<no_fm_cats; nn++) str_data_l+=strlen(focmeccats[nn]);
		if (flags.afterslip!=2) for (int nn=0; nn<all_slipmodels.NSM; nn++) str_data_l+=strlen(all_slipmodels.slipmodels[nn]);
		if (flags.afterslip!=0) for (int nn=0; nn<all_aslipmodels.NSM; nn++) str_data_l+=strlen(all_aslipmodels.slipmodels[nn]);

		str_data = malloc((str_data_l+1)* sizeof(char));
		sprintf(str_data, "%d%d%d%d%d%d%d%d%d%d", flags.err_recfault, flags.err_slipmodel,
				flags.aftershocks, (flags.aftershocks)? flags.full_field : 0,  (flags.aftershocks)? flags.aftershocks_fixedmec : 0,
				flags.afterslip, flags.OOPs, flags.err_gridpoints, use_bg_rate, use_bg_rate_file);
		sprintf(str_data, "%s%03d%03d",str_data, Nsur, Nslipmod);
		sprintf(str_data, "%s%011.5lf",str_data, tstartLL-t_firstmain);
		sprintf(str_data, "%s%s%s", str_data, crust_file, fore_template);
		sprintf(str_data, "%s%s%s%s", str_data, catname, background_rate_file);
		for (int nn=0; nn<no_fm_cats; nn++) sprintf(str_data,"%s%s",str_data, focmeccats[nn]);
		if (flags.afterslip!=2) for (int nn=0; nn<all_slipmodels.NSM; nn++) sprintf(str_data,"%s%s",str_data, all_slipmodels.slipmodels[nn]);
		if (flags.afterslip!=0) for (int nn=0; nn<all_aslipmodels.NSM; nn++) sprintf(str_data,"%s%s",str_data, all_aslipmodels.slipmodels[nn]);

		new_hash=hashlittle(str_data, str_data_l, 1);

		if (flog) {
			fprintf(flog, "\nCreated hash for this model run: %ld\n",new_hash);
			fflush(flog);
		}
	}


	//******************************************************************************//
	//  				Setup variables needed for grid search 						//
	//******************************************************************************//

	if (LLinversion && (CSEPmode || forecast)) gammas_new=dmatrix(1,Nsur,1,NgridT);
	if (LLinversion && forecast) gammas_maxLL=dmatrix(1,Nsur,1,NgridT);

	LLs=dvector(1,(1+nAsig0)*(1+nta0));
	Ldums0=dvector(1,(1+nAsig0)*(1+nta0));
	Nev=dvector(1,(1+nAsig0)*(1+nta0));
	I=dvector(1,(1+nAsig0)*(1+nta0));

	nAsig=nAsig0; nta=nta0;
	dAsig=(nAsig==0)? 0.0 : (Asig_max-Asig_min)/nAsig;	//first case to avoid 0/0 later.
	dta=(nta==0)? 0.0 : (ta_max-ta_min)/nta;

	//call these functions once over entire domain to initialize static variables in forecast_stepG2_new.
	err+=CRSLogLikelihood ((double *) 0, (double *) 0, (double *) 0, (double *)0, (double *) 0, 1, 1, DCFS, eqkfm_aft, eqkfm0res, eqkfm1, flags, Hurst,
			tevol_afterslip, crst, AllCoeff, L, max(Ntot,Nm), Nm, NgridT, focmec, fmzonelimits, NFM, &seed, cat, times2,
			fmin(tstartLL,Tstart-extra_time), tstartLL, fmax(tendLL, Tend), tw, 0.0, 0.0, r0, fixr, NULL, (double **) 0, 0, 0, 0, 0, 0, 1);

	//-----------------check if old snapshot exists, and if so read values:----------------------//

	snapshot_exists= (CSEPmode && LLinversion) ? check_if_snapshot_exists(".", &t_oldsnap, reftime, &old_hash) : 0;
	if (flog && CSEPmode) fprintf(flog, "\nSnapshot with previous LL inversion results%s found.\n", snapshot_exists? ""  : " not");
	use_snap= (snapshot_exists) ? (old_hash==new_hash) : 0;	//compare hash to verify is the same data was used in snapshot.
	if (flog && snapshot_exists) fprintf(flog, "Snapshot has %s hash as current one, and will%s be used (%ld,%ld).\n", use_snap? "same" : "different", use_snap? ""  : " not", old_hash, new_hash);

	if (use_snap) {
		if (t_oldsnap>tendLL) {
			if (verbose_level>0) printf("** Warning: old snapshot refers to a time later than tendLL (%.3lf>%.3lf)!** \n", t_oldsnap, tendLL);
			if (flog) fprintf(flog, "Warning: old snapshot refers to a time later than tendLL (%.3lf>%.3lf), and will not be used.\n", t_oldsnap, tendLL);
			use_snap=0;
		}
		else{
			err=load_oldLL(old_LLfolder, &LLs_old_matrix);
			if (err) {
				if (verbose_level>0) printf("** Warning: file containing previous LL results not found!** \n");
				if (flog) fprintf(flog, "Warning: file containing previous LL results not found.\n");
				use_snap=0;
			}
			else{
				for (int as=0; as<=nAsig; as++){
					Asig= (fixAsig)? Asig0 : Asig_min+as*dAsig;
					for (int tai=0; tai<=nta; tai++){
						ta= (fixta)? ta0 : ta_min+tai*dta;
						p+=1;
						if (LLs_old_matrix[1][p]!=Asig || LLs_old_matrix[2][p] != ta){
							if (verbose_level>0) printf("Values of (Asig, ta) from old snapshot differ from current ones - LL calculated from t=tstartLL.\n");
							if (flog) fprintf(flog, "Values of (Asig, ta) from old snapshot differ from current ones - old values will not be used.\n");
							use_snap=0;
						}
						else {
							Ldums0[p]=LLs_old_matrix[4][p];
							Nev[p]=LLs_old_matrix[5][p];
							I[p]=LLs_old_matrix[6][p];
						}
					}
				}
			}
		}
	}

	for (int p=1; p<=(1+nAsig0)*(1+nta0); p++) LLs[p]=0.0;


	//-----------------set up background rate:----------------------//

	if (flog) fprintf(flog, "\nUsing%s uniform background rate.\n", use_bg_rate? " non" : "");

	if (use_bg_rate) {
		if (use_bg_rate_file) {
			if (flog) fprintf(flog, "\nUsing background rate file %s.\n", background_rate_file);
			read_rate(crst, background_rate_file,&crst.rate0, &minmag);
			crst.r0=0;
			for (int i=1; i<=NgridT; i++) crst.r0+= crst.rate0[i]/(1.0*NgridT);
			for (int i=1; i<=NgridT; i++) crst.rate0[i]*=1.0/crst.r0;
			crst.r0*=pow(10,cat.b*(minmag-(crst.mags[1]-0.5*crst.dmags)));
			r0=crst.r0*pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));
		}
		else {
			if (use_snap) {
				sprintf(background_rate_file,"%s/%s",old_LLfolder,"backgroundrate.dat");
				if (flog) fprintf(flog, "\nUsing background rate file from previous inversion (%s).\n", background_rate_file);
				read_rate(crst, background_rate_file, &crst.rate0, NULL);
				crst.r0=0;
				for (int i=1; i<=NgridT; i++) crst.r0+= crst.rate0[i];
				for (int i=1; i<=NgridT; i++) crst.rate0[i]*=1.0/crst.r0;
				r0=crst.r0*pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));
			}
			else {
				if (flog) fprintf(flog, "\nCalculating background rate using smoothed catalog.\n");
				background_rate2(catname, &crst, reftime, 20.0, Mag_main, &(cat.Mc), &r0, 1, t_back, t_firstmain, xytoll, ztoll, smoothing, 2);
			}
		}
	}
	else {
		crst.r0=r0;
		r0=crst.r0*pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));
	}

	if (flog) fprintf(flog, "Values of background rate: \nMw>=%.2lf\t r=%.5lf\nMw>=%.2lf\t r=%.5lf\n", cat.Mc, r0, crst.mags[1]-0.5*crst.dmags, crst.r0);

	//-----------------set up LL variables:----------------------//

	if (use_snap) {
		tstartLL=t_oldsnap;
		gammas_old=dmatrix(1,Nsur,1,NgridT);
	}

	else {
		if (use_bg_rate) {
			gammas=dmatrix(0,0,1,NgridT);
			gamma_bgrate=dvector(1,NgridT);
			for (int n=1; n<=NgridT; n++) gamma_bgrate[n]=1.0/crst.rate0[n];
		}
	}

	//-----------write out snapshot for future runs of the code:------------//

	if (CSEPmode){
		fout=fopen(check_if_snapshot_filename,"w");
		times=reftime;
		times.tm_sec+=tendLL/SEC2DAY;
		mktime(&times);
		fprintf(fout, "%04d-%02d-%02dT%02d:%02d:%02dZ\n", times.tm_year+1900, times.tm_mon+1, times.tm_mday, times.tm_hour, times.tm_min, times.tm_sec);
		fprintf(fout,"%ld",new_hash);
		fclose(fout);
		if (snapshot_exists==0)	mkdir(old_LLfolder,0777);
		sprintf(fname,"%s/%s",old_LLfolder,LLsnapshot_filename);
		fout=fopen(fname,"w");
		fprintf(fout, "Asig \t ta \t r \t Ldum0 \t N \t I \t LL \n");

		if (use_bg_rate && !use_bg_rate_file){
			sprintf(fname,"%s/%s",old_LLfolder,"backgroundrate.dat");
			print_rate(fname, crst, NULL);
		}
	}
	else {
		sprintf(fname,"%s_ParamSearch.dat", outname);
		fout=fopen(fname,"w");
	}

	if (flog) fflush(flog);

	//------------------------------------------------------------------------------------------------------//
	//											 Grid search												//
	//------------------------------------------------------------------------------------------------------//

	dim=ivector(0,Nm-1);
	for (int n=0; n<Nm; n++) {
		slipmodel_combinations*=eqkfm0res[nf].parent_set_of_models->Nmod;
		dim[n]=eqkfm0res[nf].parent_set_of_models->Nmod;
		nf+=Nfaults_all[n];
	}

	maxAsig=dvector(1,slipmodel_combinations);
	maxta=dvector(1,slipmodel_combinations);
	maxr=dvector(1,slipmodel_combinations);

	if (LLinversion){
		if (verbose_level>0) printf("Performing grid search...\n");
		if (flog) fprintf(flog, "\nPerforming grid search...\nAsig \t ta \t r \t LL \n");

		//loop over all slip models:
		for (int mod=1; mod<=slipmodel_combinations; mod++){
			if (verbose_level) printf("Slip model(s) no. %d\n", mod);
			Nsm=nth_index(mod, Nm, dim);
			nf=0;
			for (int n=0; n<Nm; n++) {
				set_current_slip_model(eqkfm0res+nf,Nsm[n]);
				nf+=Nfaults_all[n];
			}
			if (flog) {
				fprintf(flog, "Using slip models:\n");
				for (int n=0; n<Nm; n++) {
					int i=0;
					for (int nn=0; nn<n; nn++) i*=dim[nn];
					i+=Nsm[n]-1;
					fprintf(flog, "\t%s\n",all_slipmodels.slipmodels[i]);
					fflush(flog);
				}
			}

			if (mod!=1 && !all_slipmodels.constant_geometry){
				free_f3tensor(AllCoeff->Coeffs_dip, 1,0,1,0,1,0);
				free_f3tensor(AllCoeff->Coeffs_st, 1,0,1,0,1,0);
				free(AllCoeff);
				setup_CoeffsDCFS(&AllCoeff, NULL, crst, eqkfm0res, eqkfm1, Nm, Ntot, Nfaults_all, which_main);
			}

			p=0;
			LLmax=-1e300;

			if (!use_snap){
				for (int p=1; p<=(nAsig+1)*(nta+1); p++)  {
					Ldums0[p]=0.0;
					Nev[p]=0.0;
					I[p]=0.0;
				}
			}

			for (int as=0; as<=nAsig; as++){
				Asig= (fixAsig)? Asig0 : Asig_min+as*dAsig;
				for (int tai=0; tai<=nta; tai++){
					ta= (fixta)? ta0 : ta_min+tai*dta;
					p+=1;
					if (use_snap){
						load_gammas(old_LLfolder, p, gammas_old, NgridT);
						gammas=gammas_old;
					}

					else {
						if (use_bg_rate) for (int i=1; i<=NgridT; i++) gammas[0][i]=(ta/Asig)*gamma_bgrate[i];
						else gammas=NULL;
					}

					err+=CRSLogLikelihood (LLs+p, Ldums0+p, Nev+p, I+p, &r, Nsur, Nslipmod, DCFS, eqkfm_aft, eqkfm0res, eqkfm1, flags, Hurst, tevol_afterslip, crst, AllCoeff,
							L, max(Nm,Ntot), Nm, NgridT, focmec, fmzonelimits, NFM, &seed, cat, times2, tstartLL, tstartLL, tendLL, tw, Asig, ta, r0, fixr, gammas,
							gammas_new, use_snap, 1, 0, 0, 0, !tai && !as);

					if (CSEPmode) write_gammas(old_LLfolder, p, gammas_new, Nsur, NgridT);

					if (LLs[p]>LLmax){
						LLmax=LLs[p];
						maxAsig[mod]=Asig;
						maxta[mod]=ta;
						maxr[mod]=r;
						if (forecast) copy_matrix(gammas_new, &gammas_maxLL, Nsur, NgridT);
					}
					fprintf(fout, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t%d\n",Asig,ta,r,Ldums0[p],Nev[p],I[p],LLs[p], mod);
					if (flog) fprintf(flog, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t%d\n", Asig, ta, r, LLs[p], mod);
					fflush(flog);
				}
			}
		}
	}

	else {
		for (int mod=1; mod<=slipmodel_combinations; mod++){
			maxta[mod]=ta0;
			maxAsig[mod]=Asig0;
			maxr[mod]=crst.r0;
			fprintf(fout, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t%d \n",maxAsig[mod],maxta[mod],0.0,0.0,0.0, mod);
		}
	}
	fclose(fout);
	if (flog) {
		fprintf(flog, "\nFinal Rate-and-State parameters:\n");
		for (int mod=1; mod<=slipmodel_combinations; mod++){
			fprintf(flog, "Slip model(s) no. %d:\t->\t", mod);
			fprintf(flog, "Asig=%.5lf \t ta=%.5lf \t r=%.5lf \n", maxAsig[mod], maxta[mod], maxr[mod]);
		}
	}

	//todo also write out LL for each slip model combination.

	//------------------------------------------------------------------------------------------------------//
	//											 Forecast													//
	//------------------------------------------------------------------------------------------------------//

	if (forecast){

		if (verbose_level>0) printf("Calculating forecast...\n");
		if (flog) fprintf(flog, "\nCalculating forecast...\n");

		crst.GRmags=assign_GRnorm(crst.mags, crst.nmags, cat.b, 1);

		if (CSEPmode && !extra_output){
			Ntts=1;
			tts=dvector(0,1);
			tts[0]=Tstart;
			tts[1]=Tend;
		}
		else {
			Ntts=ceil((Tend-Tstart)/fore_dt);
			tts=dvector(0,Ntts);
			tts[0]=Tstart;
			for (int t=1; t<=Ntts; t++) tts[t]=Tstart+fore_dt*t;
		}

		for (int mod=1; mod<=slipmodel_combinations; mod++){
				Nsm=nth_index(mod, Nm, dim);
				nf=0;
				for (int n=0; n<Nm; n++) {
					set_current_slip_model(eqkfm0res+nf,Nsm[n]);
					nf+=Nfaults_all[n];
				}
				if (flog) {
					fprintf(flog, "mod=%d. Using slip models:\n", mod);
					for (int n=0; n<Nm; n++) {
						int i=0;
						for (int nn=0; nn<n; nn++) i*=dim[nn];
						i+=Nsm[n]-1;
						fprintf(flog, "\t%s\n",all_slipmodels.slipmodels[i]);
						fflush(flog);
					}
				}

				if (!all_slipmodels.constant_geometry){
					free_f3tensor(AllCoeff->Coeffs_dip, 1,0,1,0,1,0);
					free_f3tensor(AllCoeff->Coeffs_st, 1,0,1,0,1,0);
					free(AllCoeff);
					setup_CoeffsDCFS(&AllCoeff, NULL, crst, eqkfm0res, eqkfm1, Nm, Ntot, Nfaults_all, which_main);
				}

			if (LLinversion &&  tendLL<=Tstart){
				if (flog && mod==1) fprintf(flog, "Using starting rates results from this LL inversion: ");
				gammas=gammas_maxLL;
				tstart_calc=tendLL;
				multi_gammas=1;
			}

			else{
				if (use_snap && Tstart-tendLL>t_oldsnap){
					if (flog && mod==1) fprintf(flog, "Using starting rates results an old LL inversion: ");
					p=p_found=0;
					for (int as=0; as<=nAsig && !p_found; as++){
						Asig= (fixAsig)? Asig0 : Asig_min+as*dAsig;
						for (int tai=0; tai<=nta && !p_found; tai++){
							ta= (fixta)? ta0 : ta_min+tai*dta;
							p+=1;
							if (Asig==maxAsig[mod] && ta==maxta[mod]) p_found=1;
						}
					}
					load_gammas(old_LLfolder, p, gammas_old, NgridT);
					gammas=gammas_old;
					tstart_calc=t_oldsnap;
					multi_gammas=1;
				}

				else{
					if (flog) fprintf(flog, "Using steady state starting rates: ");
					if (use_bg_rate) for (int i=1; i<=NgridT; i++) gammas[0][i]=(ta/Asig)*gamma_bgrate[i];
					else gammas=NULL;
					tstart_calc=fmin(Tstart-extra_time, t_firstmain-extra_time);
					multi_gammas=0;
				}
			}

			if (flog && mod==1) fprintf(flog, "Calculation starting time %.2lf. Forecast start time %.2lf.\n",tstart_calc, Tstart);

			if (slipmodel_combinations>1) sprintf(outnamemod,"%s%d",outname, mod);
			else sprintf(outnamemod,"%s",outname);
			if (CSEPmode && !extra_output){
				sprintf(print_forex,"%s_foremap.dat", outnamemod);
				CRSforecast(NULL, Nsur, Nslipmod, DCFS, eqkfm_aft, eqkfm0res, eqkfm1, flags, tevol_afterslip, crst, AllCoeff, L, max(Nm,Ntot), Nm, NgridT, focmec, fmzonelimits, NFM,
						&seed, cat,	times2,tstart_calc, tts, Ntts, tw, maxAsig[mod], maxta[mod], maxr[mod], gammas, multi_gammas, 1, Hurst, 0, print_forex, 0, 0, 0, 0, 0);

				if (flog) fprintf(flog, "Output file written: %s.\n",print_forex);
			}

			else {
				sprintf(print_cmb,"%s_cmbmap.dat", outnamemod);
				sprintf(print_forex,"%s_foremap.dat", outnamemod);
				sprintf(print_foret,"%s_forecast.dat", outnamemod);
				sprintf(printall_cmb,"%s_cmbmap_all.dat", outnamemod);
				sprintf(printall_forex,"%s_foremap_all.dat", outnamemod);
				sprintf(printall_foret,"%s_forecast_all", outnamemod);
				sprintf(print_LL,"%s_LLevents", outnamemod);

				CRSforecast(NULL, Nsur, Nslipmod, DCFS, eqkfm_aft, eqkfm0res, eqkfm1, flags, tevol_afterslip, crst, AllCoeff, L, max(Nm,Ntot), Nm, NgridT, focmec, fmzonelimits, NFM,
						&seed, cat, times2,tstart_calc, tts, Ntts, tw, maxAsig[mod], maxta[mod], maxr[mod], gammas, multi_gammas, 1, Hurst,
						 print_cmb, print_forex, print_foret, printall_cmb, printall_forex, printall_foret, print_LL);

				if (flog) fprintf(flog, "Output files written: %s, %s, %s, %s, %s, %s, %s.\n", print_cmb, print_forex, print_foret, printall_cmb, printall_forex, printall_foret, print_LL);

			}
		}
	}

	if (verbose_level>0) printf("Done.\n");
	if (flog) fprintf(flog, "Program completed successfully.\n");
	fclose(flog);
	return 0;

	//todo free memory.

}
