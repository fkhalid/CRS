/*
 * main.c
 *
 *  Created on: Aug 23, 2013
 *      Author: camilla
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
		//test_readZMAP_catindex();
		//background_rates();
		//test_forecast_stepG2_new();
		//test_countcolheader();
		//test_allOkada();
		// TODO: [Fahad] There should be provision for ignoring MPI
		//				  when running tests ...
		print_screen("Done!\n");
		return (0);
	}

	int err=0;

	FILE *fout, *fin, *foutfore;

	int LLinversion, forecast;	//flags indicating if parameter inversion should be done, an if forecast should be produced.
	int Nchar=120;
	char fname[Nchar],	syscopy[Nchar], outname[Nchar], outnamemod[Nchar], logfile[Nchar],
		print_LL[Nchar], print_cmb[Nchar],  print_forex[Nchar],  print_foret[Nchar],
		printall_cmb[Nchar],  printall_forex[Nchar],  printall_foret[Nchar];	//output related file names
	char infile[Nchar], fore_template[Nchar], catname[Nchar],
		background_rate_grid[Nchar], background_rate_cat[Nchar], slipmodelfile[Nchar],
		afterslipmodelfile[Nchar], modelparametersfile[Nchar], fixedmecfile[Nchar];	//input related file names
	char **focmeccats;


	//Slip model + catalog variables:
    struct slipmodels_list all_slipmodels, all_aslipmodels;	/*contain seismic and aseismic slip models respectively*/
	struct eqkfm *eqkfm_temp=0, 	//temporary structure to store earthquakes from catalog of foc. mec. (later copied in eqkfm_co)
			*eqkfm_co=0, 			//contains all slip models for seismic sources (including synthetic ones from focal mechanisms).
			*eqkfm_aft=0;			//contains all aseismic slip models.
	double Mag_main, Mc_source;	//Mag_main used to skip tw for LL calculations (todo: should not do this for forecasts).
								//Mag_main also used for declustering in case background rate; todo should put it somewhere else in input file.
								//Mc_sources is the minimum magnitude of earthquakes to be used as stress sources.

	double dt, dM, xytoll, ztoll, border;	//tolerance used to match earthquakes from catalog and from focal mechanism catalog.

	int Ntemp,		//size of eqkfm_temp
		Nco, Naf,		//number of seismic sources. size of Nfaults_co; no. of afterslips.
		*Nfaults_co=0;	//number of faults for each seismic source. size of eqkfm_co is the sum of its elements.

	double res, 	//slip model resolution
		gridresxy, 	//horizontal resolution of calculation grid
		gridresz;	//vertical resolution of calculation grid

	struct catalog cat;	//earthquake catalog
	double **focmec=0;	//matrix containing focal mechanisms
	int *fmzonelimits;	//indices of focmec corresponding to different zones
	int no_fm_cats;		//number of focal mechanisms catalogs (each defining a zone)

	//DCFS variables:
	struct pscmp *DCFS;	//contains stress fields from seismic sources (eqkfm_co)
	double dDCFS;		//min. value for which grid points should be selected for calculating stress changes from source events

	struct crust crst;	//contains description of model domain (grid) and other crust properties.
	int NgridT, 	//number of gridpoints (crst.N_allP)
		Nfocmec=0;	//number of focal mechanisms (size of focmec).


	struct Coeff_LinkList *AllCoeff;	//contains okada coefficients.

	//RateState variables:
	int fixta, fixAsig, fixr;		//flags indicating if parameters should be fixed (not inverted for).
	double Asig_min, Asig_max, Asig0;	//range of Asig for parameter search; value to be used if fixAsig=1 or LLinversion=0.
	double ta_min, ta_max, ta0;		//range of ta for parameter search; value to be used if fixta=1 or LLinversion=0.
	double r0;						//value of background rate to be used if fixr0=1 or LLinversion=0.
	int nAsig, nta;					//number of Asig, ta values to use in grid search.
	double maxAsig, maxta, maxr;	//optimal values from grid search.
	double Asig, ta, r,		// these are assigned during each grid search iteration.
			dAsig, dta;		// interval between values in grid search.

	struct flags flags;	//structure containing flags.
	int Nsur;	//no. of Monte Carlo iterations.

	struct tm reftime;	//reference time (IssueTime)
	double tstartLL, tendLL=0,	//start, end time of LL inversion.
			Tend, Tstart,	//start, end time of forecast.
			tw,				//time window to skip after each event with Mw>=Mag_main in LL calculation (due to catalog incompleteness)
			tstart_calc, 	//start of calculation time for forecast (to avoid recalculating stuff from inversion period).
			t_earliest_stress,	//time of the earliest stress change, used for forecast.
			smallest_time;	//start time of time steps for aseismic stresses, times2 (needed in rate_state_evol).
	double fore_dt;			//step between forecast output (temporal forecast).
	double *tts;			//vector containing times of forecast output (temporal forecast).
	int Ntts;				//size of tts.
	//double t_firstmain;
	double smoothing;	//smoothing used for calculating background seismicity.

	//grid search variables:
	int p=0;	//dummy variable
	int multi_gammas,	//If forecast==1, LLinversion==1 and  Tstart>=tendLL, the calculations from LL inversion will be used as starting values for forecast calculations.
						//multi_gammas is a flag indicating that this is the case.
		use_bg_rate_cat, //flag indicating that spatially variable background rate should be calculated from catalog provided.
		use_bg_rate_grid;//flag indicating that a pre-calculated spatially variable background rate should be used.

	double *LLs, *Ldums0, *Nev, *I,	//vectors containing LL inversion values for all grid search iterations.
									//LL: Log likelihood; Ldums0: first term of LL (\sum_{i} \log{\lambda(t_i, x_i)}); Nev=no. of events; I= second term of LL (integral)
									//rates in Ldums0, I should be multiplied by r(average background rate):
									//LL=
			LLmax, LL;	// maximum LL; dummy variable.

	double 	**gammas_new=NULL,		// dummy variable to store gamma values at the end of LL inversion period; overwritten at each grid search loop.
			**gammas_maxLL=NULL, 	// If forecast==1, LLinversion==1 and  Tstart>=tendLL, the calculations from LL inversion will be used as starting values for forecast calculations.
									// gammas_maxLL contains the values of gamma for all iterations, saved from LL inversion and used as starting values for forecast.
			**gammasfore=NULL;		// this will be set to gammas_maxLL if (forecast==1, LLinversion==1 and  Tstart>=tendLL), to NULL otherwise.

	//to switch between slip models.
	int slipmodel_combinations=1;	//number of slip model combinations (since each earthquake may have multiple alternative slip models).
	int *dim;	//number of slip models available for each earthquake
	int *Nsm, 	//contains the index of the slip model that should be used for each earthquake;
		nf=0;	//dummy variable to count no. of faults for each earthquake.

	//temporary variables.
	double minmag;	//minimum magnitude in background_rate_grid file.
	long seed;		//for random number generator
	int input_file_name_given=0;	//flag to check if used provided master input file.

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

	#ifdef _CRS_MPI
		// [Fahad] The file names are used in conditions in main.c for
		// 		   setting certain flags.
		MPI_Bcast(background_rate_grid,  120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(background_rate_cat,   120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(afterslipmodelfile, 	 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(fixedmecfile, 		 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(catname,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);
	#endif

//-----------------------read model parameters-------------------//

	err=read_modelparameters(modelparametersfile, &crst, reftime, &fixr, &fixAsig, &fixta, &r0, &Asig0, &ta0,
			&Asig_min, &Asig_max, &ta_min, &ta_max, &nAsig, &nta,	&tw, &fore_dt,
			&Nsur, &flags, &(cat.Mc), &Mag_main, &Mc_source, &dDCFS, &DCFS_cap,
			&dt, &dM, &xytoll, &ztoll, &border, &res, &gridresxy, &gridresz, &smoothing, &LLinversion, &forecast);

	if (err) {
		error_quit("Error reading InputModelParametersFile file %s.\n", modelparametersfile);
	}

	if (!forecast) Tend=tendLL;	//if forecast should not be produced, only needs to consider inversion time (tendLL).

//------- change flags if input files are missing:----//

	if (!focmeccats && flags.err_recfault) {
		//can not use variable receiver fault (flags.err_recfault) if catalog of focal mechanisms is not given.
		print_screen("Warning: InputCatalogFocMecFile or InputListCatalogFocMecFile not given: will not use variable receiver faults.\n");
		flags.err_recfault=0;
		flags.sources_all_iso=1;
	}

	if ((!strcmp(fixedmecfile,"")==0) && (flags.err_recfault || flags.OOPs)){
		//fixedmecfile (containing spatially variable fixed mechanisms) will not be used if variable receiver faults (focmec) or OOPS are given in parameter file flag.
		print_screen("Warning: FixedMecFile will not be used (receiver fault flag in parameter file not set to 'fixed').\n");
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
		error_quit("Errors while reading ForecastTemplate or FixedMecFile. Exiting.\n", fore_template);
	}
	NgridT=crst.N_allP;


//---------------------------------------------//
//--------------Setup afterslip----------------//
//---------------------------------------------//

	if (flags.afterslip !=0) {
		err=read_listslipmodel(afterslipmodelfile, reftime, &all_aslipmodels, res, 1);
		err+=setup_afterslip_eqkfm(all_aslipmodels, crst, &eqkfm_aft);
		if (err!=0) error_quit("Error in setting up afterslip slip model. Exiting.\n");
		Naf=all_aslipmodels.NSM;
	}
	else {
		eqkfm_aft=NULL;
		Naf=0;
	}

	//flags.splines= (flags.afterslip)? (all_aslipmodels.NSM>1): 0;	//FIXME need to do this for each afterslip
	flags.splines= (flags.afterslip)? (eqkfm_aft[0].nosnap>1): 0;

//----------------------------------------------------------//
//--------------Setup aftershocks, mainshocks --------------//
//----------------------------------------------------------//

	//read list of coseismic slip models.
	err=read_listslipmodel(slipmodelfile, reftime, &all_slipmodels, res, 0);
	if (err) error_quit("Error in reading slip model file. Exiting.\n");

	if (flags.err_recfault) {
		// catalog should be filled up to Tend (forecast end time) since it will be used to calculate LLevents for future events too
		err = setup_catalogetc(catname, focmeccats, no_fm_cats, reftime,
							   dDCFS, Mc_source, Mag_main, crst, &cat, &eqkfm_temp, &focmec, &fmzonelimits,
							   flags, &Nfocmec, &Ntemp, dt, dM,  xytoll, ztoll, border, tw,
							   tstartLL, Tend);
	}
	else {
		err = setup_catalogetc(catname, focmeccats, no_fm_cats, reftime,
							   dDCFS, Mc_source, Mag_main, crst, &cat, &eqkfm_temp,   NULL , NULL, flags,
							   NULL, &Ntemp, dt, dM,  xytoll, ztoll, border, tw,
							   tstartLL, Tend);
	}

	if (err!=0) error_quit("Error in setting up catalog. Exiting.\n");

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
		select_fm_time(focmec, &Nfocmec, Tstart);
		if (!Nfocmec) {
			print_logfile("\nNo focal mechanisms available before t=%.2lf (ForecastStartDate). Will not use multiple receiver faults.\n", Tstart);
			flags.err_recfault=0;
		}
		else {
			print_logfile("\nWill use %d receiver focal mechanisms up to t=%.2lf (ForecastStartDate)\n", Nfocmec,  Tstart);
		}
	}

//----------------------------------------------------------//
//----------------------Add slip models --------------------//
//----------------------------------------------------------//

	err=eqkfm_addslipmodels(eqkfm_temp, all_slipmodels, &eqkfm_co, Ntemp, &Nco, &Nfaults_co, dt, dM, res, crst, flags);
	if (err!=0) error_quit("Error in setting up catalog or associating events with mainshocks. Exiting.\n");



	if(LLinversion){
		print_logfile("Inversion time period: [%2.lf - %2.lf]days, ", tstartLL, tendLL);
	}

	//--------------------------------------------------------------//
	//					Setup other things							//
	//--------------------------------------------------------------//

	if (flags.err_recfault && no_fm_cats==1 && Nsur>Nfocmec) {
	    	Nsur=Nfocmec;		//it doesn't make sense to have more iterations than focal mechanisms.
			flags.sample_all=1;	//all focal mechanisms can be sampled (one per iteration).
	    }
	else flags.sample_all=0;

	if (!flags.err_recfault && flags.OOPs) flags.err_recfault=2;	//by convention, this value means OOPs when passed to calculateDCFSrandomized.
	//if (!flags.err_recfault) {
	//	if (crst.variable_fixmec){
	//		crst.nofmzones=crst.N_allP;
	//		crst.fmzone=ivector(1,crst.N_allP);
	//		for (int i=1; i<=crst.N_allP; i++) crst.fmzone[i]=i-1;
	//	}
	//	else{
	//		crst.nofmzones=1;
	//		crst.fmzone=NULL;
	//	}
	//}

	if (!crst.uniform && flags.err_gridpoints) {
			print_screen("** Warning: grid is not uniform -> grid error will not be implemented. **\n");
			print_logfile("** Warning: grid is not uniform -> grid error will not be implemented. **\n");
		flags.err_gridpoints=0;
	}

	if (!flags.err_recfault && !flags.err_gridpoints) Nsur=1;	//since there are not sources of uncertainties.

	// [Fahad] 	- We need to make sure that the number of iterations is not less
	//		   	- than the number of MPI processes; so that all processes can
	//			- contribute to the computations.
	// TODO: Rewrite the error message so that it is more helpful to the user ...
	#ifdef _CRS_MPI
		if(numProcs > Nsur) {
			if(procId == 0) {
				printf("\n Number of MPI processes: %d", numProcs);
				printf("\n Number of iterations: %d", Nsur);
			}
			error_quit("\n ** Nsur must be greater than or equal to the number of MPI processes ** \n\n");
		}
	#endif
	//--------------Setup Coefficients and DCFS struct--------------//

	#ifdef _CRS_MPI
		double coeffsStartTime, coeffsEndTime;

		MPI_Barrier(MPI_COMM_WORLD);
		coeffsStartTime = MPI_Wtime();
	#endif

	err=setup_CoeffsDCFS(&AllCoeff, &DCFS, crst, eqkfm_co, Nco, Nfaults_co, eqkfm_aft, Naf, all_aslipmodels.Nfaults);	//FIXME change 2nd last argument.

	if (err){
		error_quit("Error in setting up okada coefficients structure or associating afterslip with a mainshock.\n");
	}

	update_CoeffsDCFS(&AllCoeff, crst, eqkfm_co, Nco, Nfaults_co);

	#ifdef _CRS_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		coeffsEndTime = MPI_Wtime();

		if(procId == 0) {
			printf("\nTime - setup_CoeffsDCFS(): %f seconds\n\n", (coeffsEndTime - coeffsStartTime));
		}
	#endif

	//--------------------------------------------------------------//
	// 					Setup time steps;							//
	//--------------------------------------------------------------//

	print_screen("Setting up time steps...");

	//todo should allow this to be set from outside
	double *Cs, *ts;
	int Nfun, L=10000;	//large number needed if afterslip==1. todo check if this is the case after fixing setup_afterslip_evol.
	double *times2;

	Nfun=1;
	Cs=dvector(0,Nfun-1);
	ts=dvector(0,Nfun-1);
	Cs[0]=436.63;
	ts[0]=14.2653;

	//todo allow for general stressing history here.
	//t_earliest_stress used later to calculate tstart_calc; 1e30 ensures value is ignored (see later).
	t_earliest_stress= (Nco>0) ? eqkfm_co[0].t-1e-4 : 1e30;	//time of earliest source. fixme will change when general stressing history is allowed.

	if (flags.afterslip){

		smallest_time=fmin(t_earliest_stress, fmin(tstartLL, Tstart));	//the first time step will be before this time; this is needed in rate_state_evol.

		err=setup_afterslip_evol(smallest_time, fmax(tendLL, Tend), Cs, ts, Nfun, &eqkfm_aft,
				Naf, all_aslipmodels.Nfaults, &L, &times2, &seed);	//FIXME change input arguments
		if(err) return 1;
		print_logfile("\nSetting up time steps for calculations: %d time steps between times [%.2lf, %.2lf].\n", L, times2[0], times2[L]);

	}

	else {
		L=0;
		times2=NULL;
	}

	print_screen("done\n");


	//******************************************************************************//
	//  				Setup variables needed for grid search 						//
	//******************************************************************************//

	// if (LLinversion && forecast &&  tendLL<=Tstart) the results from LLinversion will be saved and used for forecast calculation.
	// gammas_new is overwritten for each (Asig, ta) value; the values for the optimal value are saved in gammas_maxLL.
	// TODO: these variables may be too large to be stored in memory (for large Nsur). In this case, they should not be used (and the starting rates results from LL inversion should
	// never be used in the forecast: see if condition below).
	if (LLinversion && forecast &&  tendLL<=Tstart) {
		gammas_new=dmatrix(1,Nsur,1,NgridT);
		gammas_maxLL=dmatrix(1,Nsur,1,NgridT);
	}

	LLs=dvector(1,(1+nAsig)*(1+nta));
	Ldums0=dvector(1,(1+nAsig)*(1+nta));
	Nev=dvector(1,(1+nAsig)*(1+nta));
	I=dvector(1,(1+nAsig)*(1+nta));

	dAsig=(nAsig==0)? 0.0 : (Asig_max-Asig_min)/nAsig;	//first case to avoid 0/0 later.
	dta=(nta==0)? 0.0 : (ta_max-ta_min)/nta;

	for (int p=1; p<=(1+nAsig)*(1+nta); p++) LLs[p]=0.0;

	//call these functions once over entire domain to initialize static variables in forecast_stepG2_new.
	err+=CRSLogLikelihood ((double *) 0, (double *) 0, (double *) 0, (double *)0, (double *) 0, 1, DCFS, eqkfm_aft, eqkfm_co, flags,
			crst, AllCoeff, L, Nco, Naf, NgridT, focmec, fmzonelimits, Nfocmec, &seed, cat, times2,
			fmin(tstartLL,Tstart), tstartLL, fmax(tendLL, Tend), tw, 0.0, 0.0, 0.0, r0, fixr, NULL, (double **) 0, 0, 0, 0, 1);



	//******************************************************************************//
	//  							Setup background rate 							//
	//******************************************************************************//


	print_logfile("\nUsing%s uniform background rate.\n", (use_bg_rate_cat || use_bg_rate_grid)? " non" : "");

	if (use_bg_rate_grid) {
		print_logfile("\nUsing background rate file %s.\n", background_rate_grid);
		read_rate(crst, background_rate_grid, &crst.rate0, &r0, &minmag);
		r0*=pow(10,cat.b*(minmag-cat.Mc));	//adjust rate to the catalog (since LL inversion is based on catalog).
	}
	else {
		if (use_bg_rate_cat){
			print_logfile("\nCalculating background rate using smoothed catalog from file %s.\n", background_rate_cat);
			err=background_rate(background_rate_cat, &crst, reftime, Mag_main, &minmag, &r0, &crst.rate0, xytoll, ztoll, smoothing, 2);
			if (!err) {
				r0*=pow(10,cat.b*(minmag-cat.Mc));	//adjust rate to the catalog (since LL inversion is based on catalog).
			}
			else{
				print_screen("Could not calculate background rate from smoothed catalog. will use uniform background rate.\n");
				print_logfile("Could not calculate background rate from smoothed catalog. will use uniform background rate.\n");
				use_bg_rate_cat=0;
				//settings for use_bg_rate_cat=0:
				r0*=pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));	//adjust rate to the catalog (since LL inversion is based on catalog).
				crst.rate0=NULL;	//by convention, this is equivalent to all 1s.
			}
		}
		else {
			r0*=pow(10,cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));	//adjust rate to the catalog (since LL inversion is based on catalog).
			crst.rate0=NULL;	//by convention, this is equivalent to all 1s.
		}
	}

	print_logfile("Default values of background rate: \nMw>=%.2lf\t r=%.5lf\n", cat.Mc, r0);

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

	dim=ivector(0,Nco-1);
	for (int n=0; n<Nco; n++) {
		if (eqkfm_co[nf].parent_set_of_models->Nmod) slipmodel_combinations*=eqkfm_co[nf].parent_set_of_models->Nmod;
		dim[n]=MAX(eqkfm_co[nf].parent_set_of_models->Nmod, 1);
		nf+=Nfaults_co[n];
	}

	//loop over all slip models:
	for (int mod=1; mod<=slipmodel_combinations; mod++) {
		// FIXME: [Fahad] For testing purposes only ...
		#ifdef _CRS_MPI
			dcfsStartTime = MPI_Wtime();
		#endif

		print_screen("Slip model(s) no. %d\n", mod);
		Nsm=nth_index(mod, Nco, dim);
		nf=0;
		for (int n=0; n<Nco; n++) {
			set_current_slip_model(eqkfm_co+nf,Nsm[n]);
			nf+=Nfaults_co[n];
		}

		//print information about slip models to log file:
		nf=0;
		int i=0, nn0=0; //counter: slip model names, events which have a slip model.
		print_logfile("Using slip models:\n");
		for (int n=0; n<Nco; n++) {
			//fixme this is wrong is a slip model in the list is not used.
			if (eqkfm_co[nf].parent_set_of_models->Nmod) {
				print_logfile("\t%s\n",all_slipmodels.slipmodels[i+Nsm[n]-1]);
				i+=all_slipmodels.no_slipmodels[nn0];
				nn0+=1;
			}
			else print_logfile("\t%s\n","Synthetic slip model (or isotropic field)");
			nf+=Nfaults_co[n];
		}


		if (!all_slipmodels.constant_geometry){
			if (mod!=1){
//				if (AllCoeff->Coeffs_dip) free_f3tensor(AllCoeff->Coeffs_dip, 1,0,1,0,1,0);
//				if (AllCoeff->Coeffs_st) free_f3tensor(AllCoeff->Coeffs_st, 1,0,1,0,1,0);
//				if (AllCoeff) free(AllCoeff);
				update_CoeffsDCFS(&AllCoeff, crst, eqkfm_co, Nco, Nfaults_co);
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
		maxta=ta0;
		maxAsig=Asig0;
		maxr=r0;

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

					err += CRSLogLikelihood(LLs+p, Ldums0+p, Nev+p, I+p, &r, Nsur, DCFS, eqkfm_aft,
										  	eqkfm_co, flags, crst, AllCoeff,
										  	L, Nco, Naf, NgridT, focmec, fmzonelimits, Nfocmec, &seed, cat,
										  	times2, tstartLL, tstartLL, tendLL, tw, Mag_main, Asig, ta, r0, fixr, NULL,
										  	gammas_new, 0, 0, 0, !tai && !as);

					if (!err){

						if (LLs[p]>LLmax){
							LLmax=LLs[p];
							maxAsig=Asig;
							maxta=ta;
							maxr=r;
							if (forecast &&  tendLL<=Tstart) {
								//copy gamma values into gammas_maxLL; they will be used for forecast.
								copy_matrix(gammas_new, &gammas_maxLL, Nsur, NgridT);
							}
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
				fprintf(fout, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t%d \n",maxAsig,maxta,maxr,0.0,0.0, mod);
			}
		}

		print_logfile("\nFinal values of background rate: Mw>=%.2lf\t r=%.5lf\n", cat.Mc, maxr);
		maxr*=pow(10,-cat.b*(crst.mags[1]-0.5*crst.dmags-cat.Mc));	//adjust rate to forecast magnitude range.
		print_logfile("Background rate used for forecast: Mw>=%.2lf\t r=%.5lf\n", crst.mags[1]-0.5*crst.dmags, maxr);


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
				gammasfore=gammas_maxLL;	//values from LLinversion.
				tstart_calc=tendLL;		//calculations should start from tendLL (time to which gammasfore values refer).
				multi_gammas=1;	//since gammasfore contains gamma values for each Monte Carlo iteration.

				// Set all stress fields to 0, since they will not be used anymore (but they are printed out)
				// This is done since DCFS.cmb contains the stress field from the last iteration, which will not be recalculated.
				for (int eq=0; eq<Nco; eq++){
					if (DCFS[eq].t<tstart_calc){
						for (int i=1; i<=DCFS[eq].nsel; i++){
							DCFS[eq].cmb[i]=0.0;
						}
					}
				}
			}
			else {
				if(mod==1) {
					print_logfile("Using steady state starting rates: ");
				}
				gammasfore=NULL;	// default values (will do full calculation, using uniform background rate)
				tstart_calc= fmin(Tstart, t_earliest_stress);	//start when forecast is required (Tstart) or when stress sources start t_earliest_stress.
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

			CRSforecast(&LL, Nsur, DCFS, eqkfm_aft, eqkfm_co, flags, crst, AllCoeff, L, Nco, Naf, NgridT, focmec, fmzonelimits, Nfocmec,
					&seed, cat, times2,tstart_calc, tts, Ntts, maxAsig, maxta, maxr, gammasfore, multi_gammas, 1,
					 print_cmb, print_forex, print_foret, printall_cmb, printall_forex, printall_foret, print_LL);

			print_logfile("Output files written: %s, %s, %s, %s, %s, %s, %s.\n",
								  print_cmb, print_forex, print_foret, printall_cmb, printall_forex,
								  printall_foret, print_LL);
			if(procId == 0) {
				fprintf(foutfore, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t%d\n",maxAsig, maxta, maxr, LL, mod);
			}
		}

		#ifdef _CRS_MPI
			// FIXME: [Fahad] For testing purposes only ...
			MPI_Barrier(MPI_COMM_WORLD);

			forecastEndTime = MPI_Wtime();
			forecastTotalTime += forecastEndTime - forecastStartTime;

			print_screen("\nConverting output files from Binary to ASCII ... ");

			if(procId == 0) {
				// [Fahad]: Converting file written using MPI I/O routines
				//		  : from binary to ASCII.
				char binaryToAsciiCmd_cmb[500],
					 binaryToAsciiCmd_forex[500],
					 binaryToAsciiCmd_foret_cumu[500],
					 binaryToAsciiCmd_foret_rates[500],
					 rmCmd[500], mvCmd[500];

				sprintf(binaryToAsciiCmd_cmb,
						"%s '%d/8 \"%%.6e\\t\"' -e '\"\\n\"' %s.dat > %s.txt", "hexdump -v -e",
						NgridT, printall_cmb, printall_cmb);

				sprintf(binaryToAsciiCmd_forex,
						"%s '%d/8 \"%%.6e\\t\"' -e '\"\\n\"' %s.dat > %s.txt", "hexdump -v -e",
						NgridT, printall_forex, printall_forex);

				sprintf(binaryToAsciiCmd_foret_cumu,
						"%s '%d/8 \"%%f\\t\"' -e '\"\\n\"' %s_cumu.dat > %s_cumu.txt", "hexdump -v -e",
						Ntts, printall_foret, printall_foret);

				sprintf(binaryToAsciiCmd_foret_rates,
						"%s '%d/8 \"%%f\\t\"' -e '\"\\n\"' %s_rates.dat > %s_rates.txt", "hexdump -v -e",
						Ntts, printall_foret, printall_foret);

				system(binaryToAsciiCmd_cmb);
				system(binaryToAsciiCmd_forex);
				system(binaryToAsciiCmd_foret_cumu);
				system(binaryToAsciiCmd_foret_rates);

				sprintf(mvCmd, "mv %s.txt %s.dat", printall_cmb, printall_cmb);
				system(mvCmd);

				sprintf(mvCmd, "mv %s.txt %s.dat", printall_forex, printall_forex);
				system(mvCmd);

				sprintf(mvCmd, "mv %s_cumu.txt %s_cumu.dat", printall_foret, printall_foret);
				system(mvCmd);

				sprintf(mvCmd, "mv %s_rates.txt %s_rates.dat", printall_foret, printall_foret);
				system(mvCmd);
			}

			print_screen("Done.\n");
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

	print_logfile("\nFinal Rate-and-State parameters (used for forecast):\n");
	for (int mod=1; mod<=slipmodel_combinations; mod++){
		print_logfile("Slip model(s) no. %d:\t->\t", mod);
		print_logfile("Asig=%.5lf \t ta=%.5lf \t r=%.5lf \n", maxAsig, maxta, maxr);
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
