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
#include "src/inp_out/read_crust.h"
#include "src/inp_out/read_csep_template.h"
#include "src/inp_out/read_eqkfm.h"
#include "src/inp_out/read_inputfile.h"
#include "src/seis/GR.h"
#include "src/util/hash.h"
#include "src/util/moreutil.h"
#include "src/util/nrutil.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

//todo make program work for non uniform grid (set griderr=0 in that case).
//todo check that lat=(-180, 180) are interpreted correctly throughout the program.

double DCFS_cap;
FILE *flog=NULL;
int extra_verbose, quiet;

int main (int argc, char **argv) {
	// [Fahad] Variables used by MPI related code.
	int procId = 0;		// [Fahad] Process rank
	int numProcs = 1;	// [Fahad] Total number of MPI processes
	double startTime, endTime;

	#ifdef _CRS_MPI
		omp_set_num_threads(2);

		// [Fahad] Initialize MPI and get the rank and the total number of processes.
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

		if(procId == 0) {
			printf("\n Starting CRS-MPI ... \n");
			printf("\n MPI ranks: %d", numProcs);
			printf("\n Max OMP threads per rank: %d", omp_get_max_threads());
			printf("\n Total max threads across all ranks: %d \n\n", (numProcs * omp_get_max_threads()));
		}
	#endif

	setenv("TZ", "UTC", 1);
	int run_tests=0;

	extra_verbose=0;
	quiet=0;

	if (run_tests){
		extra_verbose=1;
		//test_countcolheader();
		//test_allOkada();
		// TODO: [Fahad] There should be provision for ignoring MPI
		//				  when running tests ...
		print_screen("Done!\n");
		return (0);
	}

	int err=0;

	FILE *fout, *fin, *foutfore;

	int LLinversion, forecast;
	int Nchar=120, Nchar_long=500;
	char fname[Nchar],	infile[Nchar], fore_template[Nchar], catname[Nchar], outname[Nchar],
		syscopy[Nchar], background_rate_grid[Nchar], background_rate_cat[Nchar], slipmodelfile[Nchar], afterslipmodelfile[Nchar], modelparametersfile[Nchar],
		logfile[Nchar], print_LL[Nchar], outnamemod[Nchar];
	char print_cmb[Nchar],  print_forex[Nchar],  print_foret[Nchar],  printall_cmb[Nchar],  printall_forex[Nchar],  printall_foret[Nchar];
	char line[Nchar_long];
	char **focmeccats;
	char fixedmecfile[Nchar];

	//Slip model + catalog variables:
    struct slipmodels_list all_slipmodels, all_aslipmodels;
	struct eqkfm *eqkfm1=0, *eqkfm0res=0, *eqkfm_aft=0;	//contain respectively: aftershocks with fm, aftershocks, mainshocks, mainshocks resampled, afterslip, afterslip resampled.
	double Mag_main, Mc_source;	//Mag_main used to skip tw for LL calculations (todo: should not do this for forecasts).
								//Mag_main also used for declustering in case background rate; should put it somewhere else in input file.
	double dt, dM, xytoll, ztoll, border;

	int Nm, *Nfaults_all=0;
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

	//RateState variables
	int fixta, fixAsig, fixr;
	double Asig_min, Asig_max, Asig0, min_dAsig;
	double ta_min, ta_max, ta0, min_dta;
	double r0;
	int nAsig0, nAsig;
	int nta0, nta;
	double *maxAsig, *maxta, *maxr;
	double *rate_dum;
	double Asig, dAsig, ta, dta, r;

	struct flags flags;
	int Nsur, Nslipmod;
	int no_fm_cats;

	struct tm reftime, times;
	double tstartLL, tendLL=0, Tend, Tstart, tw, tstart_calc, tendCat, time_min_focmec; //todo read tendLL from file.
	double fore_dt;
	double *tts;
	int Ntts;
	double t_firstmain;
	double smoothing;	//for calculating background seismicity.

	//grid search variables:
	int p=0, p_found;
	int multi_gammas, use_bg_rate_cat, use_bg_rate_grid;
	double *LLs, *Ldums0, *Nev, *I, LLmax, LL;
	double *gamma_bgrate=NULL;
	double 	*gammas=NULL;
	double 	**gammas_new=NULL, \
			**gammas_old=NULL, \
			**gammas_maxLL=NULL, \
			**gammasfore=NULL;

	//to switch between slip models.
	int refresh;
	int slipmodel_combinations=1;
	int *dim;
	int *Nsm, nf=0;

	//temporary variables.
	double minmag;
	long seed;
	int nc;
	int j0, N;
	int input_file_name_given=0;

	// FIXME: [Fahad] For testing purposes only ...
	#ifdef _CRS_MPI
		startTime = MPI_Wtime();
	#endif


	//-----------------------Parse input from command line --------------------//

	if(procId == 0) {
		if (argc==1) {
			print_screen("Usage: %s input_file [options]\nOptions:\n\t --verbose (-v) : outputs extra messages to logfile and screen;\n\t --quiet (-q) : no output.\n\n", argv[0]);
			return 1;
		}

		for (int i=1; i<argc; i++){
			if (!strncmp(argv[i],"-",1)){
				if (!(strcmp(argv[i],"--verbose")) || !(strcmp(argv[i],"-v"))) {
					if (quiet) {
						error_quit("Error: options --verbose (-v) and --quiet (-q) are incompatible!\n\n");
					}
					else extra_verbose=1;
				}
				else {
					if (!(strcmp(argv[i],"--quiet")) || !(strcmp(argv[i],"-q"))) {
						if (extra_verbose) {
							error_quit("Error: options --verbose (-v) and --quiet (-q) are incompatible!\n\n");
						}
						else quiet=1;
					}
					else {
						print_screen("Invalid option %s\n\n", argv[i]);
						return 1;
					}
				}
			}
			else{
				input_file_name_given=1;
				print_screen("Input file: %s\n", argv[1]);
				sscanf(argv[i],"%s", infile);
			}
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&input_file_name_given, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (!input_file_name_given){
		print_screen("Usage: %s input_file [options]\nOptions:\n\t --verbose (-v) : outputs extra messages to logfile and screen;\n\t --quiet (-q) : no output.\n\n", argv[0]);
		return 1;
	}

	//-----------------------read input file -------------------//

	err=read_inputfile(infile, outname, fore_template, catname, &focmeccats, background_rate_grid, background_rate_cat,
			fixedmecfile, slipmodelfile, afterslipmodelfile, modelparametersfile, logfile, &reftime, &Tstart, &Tend, &tstartLL, &seed,
			&no_fm_cats);

	if (err) {
		error_quit("Error reading input file %s.\n", infile);
	}

	//todo [askFahad]: why does this have to be repeated? (also in read_inputfile).
						//(and why is focmecfile not repeated?)
	#ifdef _CRS_MPI
		// [Fahad] The file names are used in conditions in main.c for
		// 		   setting certain flags.
		MPI_Bcast(background_rate_grid,  120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(background_rate_cat,   120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(afterslipmodelfile, 	 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(fixedmecfile, 		 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(catname,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);	//todo [askFahad]: do we need to broadcast this?
	#endif

//-----------------------read model parameters-------------------//

	err=read_modelparameters(modelparametersfile, &crst, reftime, &fixr, &fixAsig, &fixta, &r0, &Asig0, &ta0,
			&Asig_min, &Asig_max, &ta_min, &ta_max, &nAsig0, &nta0,	&tw, &fore_dt,
			&Nsur, &flags, &(cat.Mc), &Mag_main, &Mc_source, &dDCFS, &DCFS_cap,
			&dt, &dM, &xytoll, &ztoll, &border, &res, &gridresxy, &gridresz, &smoothing, &LLinversion, &forecast);

	if (err) {
		error_quit("Error reading InputModelParametersFile file %s.\n", modelparametersfile);
	}
	
	//future events (cat.t[i]>0) are only needed to calculate LL for forecast, if forecast is produced.
	//todo remove if decides not to print out LLevents.
	tendCat= (forecast)? Tend : 0;

//------- change flags if input files are missing:----//

	if (!focmeccats && flags.err_recfault) {
		print_screen("Warning: InputCatalogFocMecFile or InputListCatalogFocMecFile not given: will not use variable receiver faults.\n");
		flags.err_recfault=0;
		flags.sources_all_iso=1;
	}
	if ((strcmp(afterslipmodelfile,"")==0)) {
		print_screen("InputListAfterslipModels not given: will not use afterslip.\n");
		flags.afterslip=0;
	}
	else{
		flags.afterslip=1;
	}

	if (strcmp(background_rate_grid,"")==0)	use_bg_rate_grid=0;
	else use_bg_rate_grid=1;

	if (strcmp(background_rate_cat,"")==0)	use_bg_rate_cat=0;
	else use_bg_rate_cat=1;

//----------- Copy input and parameters file to log file -------//

	if(procId == 0) {
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
		}
		else flog=NULL;

		sprintf(syscopy,"cp %s %s_inputfile.txt",infile, outname);
		system(syscopy);
		sprintf(syscopy,"cp %s %s_parameters.txt", modelparametersfile, outname);
		system(syscopy);
	}

//----------- Read grid (crst) information from templates -------//

	err=read_crust(fore_template, fixedmecfile , &crst, gridresxy, gridresz, flags.err_recfault);

	if (err) {
		error_quit("Errors while reading template file %s. Exiting.", fore_template);
	}
	NgridT=crst.N_allP;


//---------------------------------------------//
//--------------Setup afterslip----------------//
//---------------------------------------------//

	if (flags.afterslip !=0) {
		read_listslipmodel(afterslipmodelfile, reftime, &all_aslipmodels, res, 1);
		err=setup_afterslip_eqkfm(all_aslipmodels, crst, &eqkfm_aft);
		if (err!=0) error_quit("Error in setting up afterslip slip model - exiting.");
	}
	else eqkfm_aft=NULL;

	flags.splines= (flags.afterslip)? (all_aslipmodels.NSM>1): 0;

//----------------------------------------------------------//
//--------------Setup aftershocks, mainshocks --------------//
//----------------------------------------------------------//

	//read list of coseismic slip models.
	err=read_listslipmodel(slipmodelfile, reftime, &all_slipmodels, res, 0);
	if (err) error_quit("Error in reading slip model file. Exiting.\n");

	if (flags.err_recfault) {
		err = setup_catalogetc(catname, focmeccats, no_fm_cats, reftime,
							   dDCFS, Mc_source, Mag_main, crst, &cat, &eqkfm1, &focmec, &fmzonelimits,
							   flags, &NFM, &Ntot, dt, dM,  xytoll, ztoll, border, tw,
							   tstartLL, tendCat);
	}
	else {
		err = setup_catalogetc(catname, focmeccats, no_fm_cats, reftime,
							   dDCFS, Mc_source, Mag_main, crst, &cat, &eqkfm1,   NULL , NULL, flags,
							   NULL, &Ntot, dt, dM,  xytoll, ztoll, border, tw,
							   tstartLL, tendCat);
	}

	if (err!=0) error_quit("**Error in setting up catalog. Exiting. **");

	if (flags.err_recfault && (no_fm_cats!=crst.nofmzones)){
		if (crst.nofmzones>no_fm_cats){
			print_screen("**Error: not enough catalogs of focal mechanisms given! (%d given, at least %d required)**\n", no_fm_cats, crst.nofmzones);
			print_logfile("**Error: not enough catalogs of focal mechanisms given! (%d given, at least %d required)**\n", no_fm_cats, crst.nofmzones);
			crst.nofmzones=1;
		}
		else {
			print_screen("**Warning: some catalogs of focal mechanisms not used! (%d given, %d used)**\n", no_fm_cats, crst.nofmzones);
			print_logfile("**Warning: some catalogs of focal mechanisms not used! (%d given, %d used)**\n", no_fm_cats, crst.nofmzones);
		}
	}

//----------------------------------------------------------//
//-----------------Setup LL inversion period ---------------//
//----------------------------------------------------------//

	if (flags.err_recfault){
		//only use focal mechanisms before start of LL period (Tstart).
		select_fm_time(focmec, &NFM, Tstart);
		if (!NFM) {
			print_logfile("\nNo focal mechanisms available before t=%.2lf (ForecastStartDate). Will not use multiple receiver faults.\n", Tstart);
			flags.err_recfault=0;
		}
		else {
			print_logfile("\nWill use %d receiver focal mechanisms up to t=%.2lf (ForecastStartDate)\n", NFM,  Tstart);
		}
	}

//----------------------------------------------------------//
//----------------------Add slip models --------------------//
//----------------------------------------------------------//

	err=eqkfm_addslipmodels(eqkfm1, all_slipmodels, &eqkfm0res, Ntot, &Nm, &Nfaults_all, dt, dM, res, crst, flags);
	if (err!=0) error_quit("**Error in setting up catalog or associating events with mainshocks. Exiting. **\n");

	if(LLinversion){
		print_logfile("Inversion time period: [%2.lf - %2.lf]days, ", tstartLL, tendLL);
	}

	//--------------Setup Coefficients and DCFS struct--------------//

	#ifdef _CRS_MPI
		double coeffsStartTime, coeffsEndTime;

		MPI_Barrier(MPI_COMM_WORLD);
		coeffsStartTime = MPI_Wtime();
	#endif

	if (flags.afterslip){
		err=setup_CoeffsDCFS(&AllCoeff, &DCFS, crst, eqkfm0res, Nm, Nfaults_all, all_aslipmodels.tmain[-1], 1);
	}
	else{
		err=setup_CoeffsDCFS(&AllCoeff, &DCFS, crst, eqkfm0res, Nm, Nfaults_all, 0.0, 0);
	}
	if (err){
		error_quit("Error in setting up okada coefficients structure or associating afterslip with a mainshock.\n");
	}

	update_CoeffsDCFS(&AllCoeff, crst, eqkfm0res, Nm, Nfaults_all);

	#ifdef _CRS_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		coeffsEndTime = MPI_Wtime();

		if(procId == 0) {
			printf("\nTime - setup_CoeffsDCFS(): %f seconds\n\n", (coeffsEndTime - coeffsStartTime));
		}
	#endif

	//--------------------------------------------------------------//
	//					Setup other things							//
	//--------------------------------------------------------------//

	if (!flags.err_recfault && !flags.err_gridpoints) Nsur=1;	//since there are not sources of uncertainties.
	if (flags.err_recfault && no_fm_cats==1 && Nsur>NFM) {
	    	Nsur=NFM;
			flags.sample_all=1;
	    }
	else flags.sample_all=0;

	if (!flags.err_recfault && flags.OOPs) flags.err_recfault=2;	//by convention, this value means OOPs when passed to calculateDCFSrandomized.
	if (!flags.err_recfault) {
		if (crst.variable_fixmec){
			crst.nofmzones=crst.N_allP;
			crst.fmzone=ivector(1,crst.N_allP);
			for (int i=1; i<=crst.N_allP; i++) crst.fmzone[i]=i-1;
		}
		else{
			crst.nofmzones=1;
			crst.fmzone=NULL;
		}
	}

	if (!crst.uniform && flags.err_gridpoints) {
			print_screen("** Warning: grid is not uniform -> grid error will not be implemented. **\n");
			print_logfile("** Warning: grid is not uniform -> grid error will not be implemented. **\n");
		flags.err_gridpoints=0;
	}

	//--------------------------------------------------------------//
	// 					Setup time steps;							//
	//--------------------------------------------------------------//

	print_screen("Setting up time steps...");

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

	print_screen("done\n");
	print_logfile("\nSetting up time steps for calculations: %d time steps between times [%.2lf, %.2lf].\n", L, times2[0], times2[L]);


	//******************************************************************************//
	//  				Setup variables needed for grid search 						//
	//******************************************************************************//

	if (LLinversion && forecast) {
		gammas_new=dmatrix(1,Nsur,1,NgridT);
		gammas_maxLL=dmatrix(1,Nsur,1,NgridT);
	}

	LLs=dvector(1,(1+nAsig0)*(1+nta0));
	Ldums0=dvector(1,(1+nAsig0)*(1+nta0));
	Nev=dvector(1,(1+nAsig0)*(1+nta0));
	I=dvector(1,(1+nAsig0)*(1+nta0));

	nAsig=nAsig0; nta=nta0;
	dAsig=(nAsig==0)? 0.0 : (Asig_max-Asig_min)/nAsig;	//first case to avoid 0/0 later.
	dta=(nta==0)? 0.0 : (ta_max-ta_min)/nta;

	//call these functions once over entire domain to initialize static variables in forecast_stepG2_new.
	err+=CRSLogLikelihood ((double *) 0, (double *) 0, (double *) 0, (double *)0, (double *) 0, 1, DCFS, eqkfm_aft, eqkfm0res, flags,
			tevol_afterslip, crst, AllCoeff, L, Nm, NgridT, focmec, fmzonelimits, NFM, &seed, cat, times2,
			fmin(tstartLL,Tstart), tstartLL, fmax(tendLL, Tend), tw, 0.0, 0.0, 0.0, r0, fixr, NULL, (double **) 0, 0, 0, 0, 1);

	for (int p=1; p<=(1+nAsig0)*(1+nta0); p++) LLs[p]=0.0;


	//-----------------set up background rate:----------------------//

	print_logfile("\nUsing%s uniform background rate.\n", (use_bg_rate_cat || use_bg_rate_grid)? " non" : "");

	if (use_bg_rate_grid) {
		print_logfile("\nUsing background rate file %s.\n", background_rate_grid);
		read_rate(crst, background_rate_grid,&crst.rate0, &minmag);
		crst.r0=0;
		for (int i=1; i<=NgridT; i++) crst.r0+= crst.rate0[i];
		for (int i=1; i<=NgridT; i++) crst.rate0[i]*=crst.N_allP/crst.r0;
		crst.r0*=pow(10,cat.b*(minmag-(crst.mags[1]-0.5*crst.dmags)));
		r0=crst.r0*pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));
	}
	else {
		if (use_bg_rate_cat){
			print_logfile("\nCalculating background rate using smoothed catalog from file %s.\n", background_rate_cat);
			err=background_rate2(background_rate_cat, &crst, reftime, 20.0, Mag_main, &(cat.Mc), &r0, 1, xytoll, ztoll, smoothing, 2);
			crst.r0=r0*pow(10,cat.b*(cat.Mc-crst.mags[1]+0.5*crst.dmags));
			if (err){
				print_screen("Could not calculate background rate from smoothed catalog. will use uniform background rate.\n");
				print_logfile("Could not calculate background rate from smoothed catalog. will use uniform background rate.\n");
				crst.r0=r0;
				r0=crst.r0*pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));
				use_bg_rate_cat=0;
				crst.rate0=NULL;	//by convention, this is equivalent to all 1s.
			}
		}
		else {
			crst.r0=r0;
			r0=crst.r0*pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));
			crst.rate0=NULL;	//by convention, this is equivalent to all 1s.
		}
	}

	print_logfile("Values of background rate: \nMw>=%.2lf\t r=%.5lf\nMw>=%.2lf\t r=%.5lf\n", cat.Mc, r0, crst.mags[1]-0.5*crst.dmags, crst.r0);

	//-----------------set up LL variables:----------------------//

	if (use_bg_rate_cat || use_bg_rate_grid) {
		gammas=dvector(1,NgridT);
		gamma_bgrate=dvector(1,NgridT);
		for (int n=1; n<=NgridT; n++) gamma_bgrate[n]=1.0/crst.rate0[n];
	}


	//-----------write out summary of grid search:------------//

	if(procId == 0) {
		sprintf(fname,"%s_ParamSearch.dat", outname);
		fout=fopen(fname,"w");
		sprintf(fname,"%s_LogLikelihood.dat", outname);
		foutfore=fopen(fname,"w");
	}

	//-----------Setup variable needed for forecast:------------//

	if (forecast){
		crst.GRmags=assign_GRnorm(crst.mags, crst.nmags, cat.b, 1);
		Ntts=ceil((Tend-Tstart)/fore_dt);
		tts=dvector(0,Ntts);
		tts[0]=Tstart;
		for (int t=1; t<=Ntts; t++) tts[t]=Tstart+fore_dt*t;
	}


	//------------------------------------------------------------------------------------------------------//
	//								  Grid search and forecast												//
	//------------------------------------------------------------------------------------------------------//

	#ifdef _CRS_MPI
		// [Fahad] Make sure all processes are in sync at this point.
		MPI_Barrier(MPI_COMM_WORLD);

		// Timing for benchmarking
		endTime = MPI_Wtime();
		if(procId == 0) {
			printf("\nTime - I/O + broadcast: %f seconds\n\n", (endTime - startTime));
		}

		startTime = MPI_Wtime();

		double dcfsStartTime, dcfsEndTime, dcfsTotalTime = 0.0;
		double gridStartTime, gridEndTime, gridTotalTime = 0.0;
		double forecastStartTime, forecastEndTime, forecastTotalTime = 0.0;
	#endif

	dim=ivector(0,Nm-1);
	for (int n=0; n<Nm; n++) {
		if (eqkfm0res[nf].parent_set_of_models->Nmod) slipmodel_combinations*=eqkfm0res[nf].parent_set_of_models->Nmod;
		dim[n]=MAX(eqkfm0res[nf].parent_set_of_models->Nmod, 1);
		nf+=Nfaults_all[n];
	}

	maxAsig=dvector(1,slipmodel_combinations);
	maxta=dvector(1,slipmodel_combinations);
	maxr=dvector(1,slipmodel_combinations);

	//loop over all slip models:
	for (int mod=1; mod<=slipmodel_combinations; mod++) {
		// FIXME: [Fahad] For testing purposes only ...
		#ifdef _CRS_MPI
			dcfsStartTime = MPI_Wtime();
		#endif

		print_screen("Slip model(s) no. %d\n", mod);
		Nsm=nth_index(mod, Nm, dim);
		nf=0;
		for (int n=0; n<Nm; n++) {
			set_current_slip_model(eqkfm0res+nf,Nsm[n]);
			nf+=Nfaults_all[n];
		}

		//print information about slip models to log file:
		nf=0;
		int i=0, nn0=0; //counter: slip model names, events which have a slip model.
		print_logfile("Using slip models:\n");
		for (int n=0; n<Nm; n++) {
			//fixme this is wrong is a slip model in the list is not used.
			if (eqkfm0res[nf].parent_set_of_models->Nmod) {
				print_logfile("\t%s\n",all_slipmodels.slipmodels[i+Nsm[n]-1]);
				i+=all_slipmodels.no_slipmodels[nn0];
				nn0+=1;
			}
			else print_logfile("\t%s\n","Synthetic slip model (or isotropic field)");
			nf+=Nfaults_all[n];
		}


		if (!all_slipmodels.constant_geometry){
			if (mod!=1){
//				if (AllCoeff->Coeffs_dip) free_f3tensor(AllCoeff->Coeffs_dip, 1,0,1,0,1,0);
//				if (AllCoeff->Coeffs_st) free_f3tensor(AllCoeff->Coeffs_st, 1,0,1,0,1,0);
//				if (AllCoeff) free(AllCoeff);
				update_CoeffsDCFS(&AllCoeff, crst, eqkfm0res, Nm, Nfaults_all);
			}
		}

		// FIXME: [Fahad] For testing purposes only ...
		#ifdef _CRS_MPI
			MPI_Barrier(MPI_COMM_WORLD);

			dcfsEndTime = MPI_Wtime();
			dcfsTotalTime += dcfsEndTime - dcfsStartTime;
		#endif

		//------------------------------------------//
		//				Grid Search					//
		//------------------------------------------//

		// FIXME: [Fahad] For testing purposes only ...
		#ifdef _CRS_MPI
			gridStartTime = MPI_Wtime();
		#endif

		//set default values:
		maxta[mod]=ta0;
		maxAsig[mod]=Asig0;
		maxr[mod]=crst.r0;

		if (LLinversion) {
			print_screen("Performing grid search...\n");
			print_logfile("\nPerforming grid search...\nAsig \t ta \t r \t LL \n");

			p=0;
			LLmax=-1e300;

			for (int p=1; p<=(nAsig+1)*(nta+1); p++)  {
				Ldums0[p]=0.0;
				Nev[p]=0.0;
				I[p]=0.0;
			}

			for(int as=0; as<=nAsig; as++) {
				Asig= (fixAsig)? Asig0 : Asig_min+as*dAsig;
				for(int tai=0; tai<=nta; tai++) {
					err=0;
					ta= (fixta)? ta0 : ta_min+tai*dta;
					p+=1;
					gammas=NULL;	//fixme check if this should be allowed to be set equal to gammas_bg_rate?

					err += CRSLogLikelihood(LLs+p, Ldums0+p, Nev+p, I+p, &r, Nsur, DCFS, eqkfm_aft,
										  	eqkfm0res, flags, tevol_afterslip, crst, AllCoeff,
										  	L, Nm, NgridT, focmec, fmzonelimits, NFM, &seed, cat,
										  	times2, tstartLL, tstartLL, tendLL, tw, Mag_main, Asig, ta, r0, fixr, gammas,
										  	gammas_new, 0, 0, 0, !tai && !as);

					if (!err){

						if (LLs[p]>LLmax){
							LLmax=LLs[p];
							maxAsig[mod]=Asig;
							maxta[mod]=ta;
							maxr[mod]=r;
							if (forecast) copy_matrix(gammas_new, &gammas_maxLL, Nsur, NgridT);
						}
						if(procId == 0) {
							fprintf(fout, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t%d\n",
									Asig,ta,r,Ldums0[p],Nev[p],I[p],LLs[p], mod);
							if (flog) fprintf(flog, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t%d\n", Asig, ta, r, LLs[p], mod);
							fflush(flog);
						}
					}
					else{
						if(procId == 0) {
							fprintf(fout, "%.5lf \t %.5lf \t %.5lf \t NaN \t NaN \t NaN \t NaN \t%d\n",Asig,ta,r,mod);
							if (flog) fprintf(flog, "%.5lf \t %.5lf \t %.5lf \t NaN \t%d\n", Asig, ta, r, mod);
							fflush(flog);
						}
					}
				}
			}
		}
		else {
			if(procId == 0) {
				fprintf(fout, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t%d \n",maxAsig[mod],maxta[mod],maxr[mod],0.0,0.0, mod);
			}
		}

		//------------------------------------------//
		//			 	Forecast					//
		//------------------------------------------//

		#ifdef _CRS_MPI
			// [Fahad] Make sure all processes are in synch at this point.
			MPI_Barrier(MPI_COMM_WORLD);

			gridEndTime = MPI_Wtime();
			gridTotalTime += gridEndTime - gridStartTime;

			forecastStartTime = MPI_Wtime();
		#endif

		if (forecast) {
			print_screen("Calculating forecast...\n");
			print_logfile("\nCalculating forecast...\n");

			if (LLinversion &&  tendLL<=Tstart){
				if(mod==1) {
					print_logfile("Using starting rates results from LL inversion: ");
				}
				gammasfore=gammas_maxLL;
				tstart_calc=tendLL;
				multi_gammas=1;
			}
			else {
				print_logfile("Using steady state starting rates: ");
				if (use_bg_rate_cat || use_bg_rate_grid) for (int i=1; i<=NgridT; i++) gammasfore[0][i]=(ta/Asig)*gamma_bgrate[i];
				else gammasfore=NULL;
				tstart_calc=fmin(Tstart, t_firstmain);
				multi_gammas=0;
			}

			if(mod==1) {
				print_logfile("Calculation starting time %.2lf. Forecast start time %.2lf.\n",tstart_calc, Tstart);
			}

			if (slipmodel_combinations>1) sprintf(outnamemod,"%s%d",outname, mod);
			else sprintf(outnamemod,"%s",outname);

			sprintf(print_cmb,"%s_cmbmap", outnamemod);
			sprintf(print_forex,"%s_foremap", outnamemod);
			sprintf(print_foret,"%s_forecast", outnamemod);
			sprintf(printall_cmb,"%s_cmbmap_all", outnamemod);
			sprintf(printall_forex,"%s_foremap_all", outnamemod);
			sprintf(printall_foret,"%s_forecast_all", outnamemod);
			sprintf(print_LL,"%s_LLevents", outnamemod);

			CRSforecast(&LL, Nsur, DCFS, eqkfm_aft, eqkfm0res, flags, tevol_afterslip, crst, AllCoeff, L, Nm, NgridT, focmec, fmzonelimits, NFM,
					&seed, cat, times2,tstart_calc, tts, Ntts, maxAsig[mod], maxta[mod], maxr[mod], gammasfore, multi_gammas, 1,
					 print_cmb, print_forex, print_foret, printall_cmb, printall_forex, printall_foret, print_LL);

			print_logfile("Output files written: %s, %s, %s, %s, %s, %s, %s.\n",
								  print_cmb, print_forex, print_foret, printall_cmb, printall_forex,
								  printall_foret, print_LL);
			if(procId == 0) {
				fprintf(foutfore, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t%d\n",maxAsig[mod], maxta[mod], maxr[mod], LL, mod);
			}
		}

		#ifdef _CRS_MPI
			// FIXME: [Fahad] For testing purposes only ...
			MPI_Barrier(MPI_COMM_WORLD);

			forecastEndTime = MPI_Wtime();
			forecastTotalTime += forecastEndTime - forecastStartTime;
		#endif
	}

	// FIXME: [Fahad] For testing purposes only ...
	#ifdef _CRS_MPI
		MPI_Barrier(MPI_COMM_WORLD);

		endTime = MPI_Wtime();

		if(procId == 0) {
			printf("\n\nTime - DCFS: %f seconds", dcfsTotalTime);
			printf("\nTime - Grid Search: %f seconds", gridTotalTime);
			printf("\nTime - Forecast: %f seconds", forecastTotalTime);
			printf("\nTime - DCFS + Grid Search + Forecast: %f seconds\n\n", (endTime - startTime));
		}
	#endif

	if(procId == 0) {
		fclose(fout);
		fclose(foutfore);
	}

	print_logfile("\nFinal Rate-and-State parameters:\n");
	for (int mod=1; mod<=slipmodel_combinations; mod++){
		print_logfile("Slip model(s) no. %d:\t->\t", mod);
		print_logfile("Asig=%.5lf \t ta=%.5lf \t r=%.5lf \n", maxAsig[mod], maxta[mod], maxr[mod]);
	}

	print_screen("Done.\n");
	print_logfile("Program completed successfully.\n");
	if (flog){
		sprintf(syscopy,"date >> %s", logfile);
		system(syscopy);
		fclose(flog);
	}


	#ifdef _CRS_MPI
		MPI_Finalize();
	#endif

	return 0;

	//todo free memory.

}
