/*
 * read_param.c
 *
 *  Created on: Dec 20, 2013
 *      Author: camcat
 */

#include <stdio.h>
#include <time.h>
#include "../defines.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int read_modelparameters(char *modelparametersfile, struct crust *crst, struct tm reftime,
						int *fixr, int *fixAsig, int *fixta, double *r0, double *Asig0,
						double *ta0, double *Asig_min, double *Asig_max, double *ta_min,
						double *ta_max, int *nAsig0, int *nta0, double *tw, double *fore_dt,
						int *Nsur, struct flags *flags,
						double *Mc, double *Mag_main, double *Mc_source,
						double *dCFS, double *DCFS_cap, double *dt, double *dM,
						double *xytoll, double *ztoll, double *border, double *res,
						double *gridresxy, double *gridresz, double *smoothing,
						int *LLinversion, int *forecast) {
	// [Fahad] Variables used for MPI.

	int procId = 0;
	int fileError = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	FILE * fin;
	char comment[]="#", comm=comment[0];
	int Nchar_long=500;
	char line[Nchar_long];
	char sourcemode_fm[10], sourcemode_nofm[10];
	char recfault[10];
	struct tm times;
	char regstress_mode[120];
	double s[3];	//regional stress field description;
	double st[3];	//regional stress field description;
	double di[3];	//regional stress field description;

	// [Fahad] If there is a file error, only root will know about it.
	//		   So it is important that the error is broadcast to all
	//		   processes, and they all return from the function with
	//		   the same error code.
	if(procId == 0) {
		fin = fopen(modelparametersfile,"r");
		if(fin == NULL) {
			print_screen("Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);
			print_logfile("Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	//initialize crust structure:
	init_crst(crst);

	if(procId == 0) {
		sprintf(comment,"#");

		//-------------Coulomb stress parameters------------------//
		comm=comment[0];
		line[0]=comm;
		while (line[0]==comm)fgets(line,Nchar_long,fin);
		sscanf(line,"%lf %lf", Mc_source, border);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%s %s", sourcemode_fm, sourcemode_nofm);

		if (!strcmp(sourcemode_fm, "iso")){
			(*flags).sources_all_iso=1;
			(*flags).sources_without_focmec=1;
		}
		else {
			if (!strcmp(sourcemode_fm, "fm")){
				(*flags).sources_all_iso=0;
				if (!strcmp(sourcemode_nofm, "no")){
					(*flags).sources_without_focmec=0;
				}
				else{
					if (!strcmp(sourcemode_nofm, "iso")){
						(*flags).sources_without_focmec=1;
					}
					else{
						if (!strcmp(sourcemode_nofm, "fix")){
							(*flags).sources_without_focmec=2;
						}
						else{
							print_screen("Illegal value for source_mode_nofocmec (should be one of: no/iso/fix).\n", modelparametersfile);
							print_logfile("Illegal value for source_mode_nofocmec (should be one of: no/iso/fix).\n", modelparametersfile);
							fileError=1;
						}
					}
				}
			}
			else {
				print_screen("Illegal value for source_mode_focmec (should be one of: iso/fm).\n", modelparametersfile);
				print_logfile("Illegal value for source_mode_focmec (should be one of: iso/fm).\n", modelparametersfile);
				fileError=1;
			}
		}

		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", res);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", dCFS, DCFS_cap);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", gridresxy, gridresz);

		//-------------Parameter inversion------------------//

		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%d", LLinversion);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf", fixr, r0);
		//the following 2 lines (Asig, ta values) can have 2 alternative forms:
		//1 Asig0
		//or:
		//0 Asig1 Asig2
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf %lf  %d", fixAsig, Asig_min, Asig_max, nAsig0);
		if (*fixAsig && *LLinversion){
			*Asig0=*Asig_min;
			*Asig_min=*Asig_max=0.0;
			*nAsig0=0;
		}
		else{
			*Asig0=0.0;
		}
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf %lf %d", fixta, ta_min, ta_max, nta0);
		if (*fixta && *LLinversion){
			*ta0=*ta_min;
			*ta_min=*ta_max=0.0;
			*nta0=0;
		}
		else{
			*ta0=0.0;
		}
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", Mc);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", tw, Mag_main);

		//-------------Treatment of uncertainties------------------//
		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%s", recfault);
		if (!strcmp(recfault, "fixed")){
			(*flags).err_recfault=0;
			(*flags).OOPs=0;
		}
		else{
			if (!strcmp(recfault,"oops")){
				(*flags).err_recfault=0;
				(*flags).OOPs=1;
			}
			else {
				if (!strcmp(recfault,"focmec")){
					(*flags).err_recfault=1;
					(*flags).OOPs=0;
				}
				else{
					print_screen("Illegal value for choice of receiver fault in file %s (should be one of: oops/fixed/focmec).\n", modelparametersfile);
					print_logfile("Illegal value for choice of receiver fault in file %s (should be one of: oops/fixed/focmec).\n", modelparametersfile);
					fileError=1;
				}
			}

		}
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", &((*flags).err_gridpoints));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", Nsur);


		//-------------Other parameters------------------//

		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%d %lf", forecast, fore_dt);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", dt);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", dM);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", xytoll);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", ztoll);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", smoothing);


		//-------------Parameters Describing crustal properties------------------//

		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%lf %lf", &((*crst).lambda), &((*crst).mu));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", &((*crst).fric), &((*crst).skepton));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf %lf", (*crst).str0, (*crst).dip0, (*crst).rake0);	//todo explain somewhere the 2 ways this variables can be used.
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%s", regstress_mode);

		if (!(strcmp(regstress_mode,"oops"))) {
			fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
			sscanf(line,"%lf %lf %lf", s, s+1, s+2);
		}
		else {
			if (!(strcmp(regstress_mode,"paxis"))) {
				fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%lf %lf %lf", s, st, di);
				fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%lf %lf %lf", s+1, st+1, di+1);
				fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%lf %lf %lf", s+2, st+2, di+2);
			}

			else {
				print_logfile("Invalid value for mode of regional stress field: should be 'oops' or 'paxis'. Exit.\n");
				print_screen("Invalid value for mode of regional stress field: should be 'oops' or 'paxis'. Exit.\n");
				fileError=1;
			}
		}

		fclose(fin);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	// TODO: [Fahad] Look into how much of the following code can be moved to
	//		 a separate function ...




	#ifdef _CRS_MPI
		// Copy scalars to the BCast_Model_Parameters struct
		struct BCast_Model_Parameters modelParams;

		if(procId == 0) {
			modelParams.fixr         =  *fixr;
			modelParams.fixAsig      =  *fixAsig;
			modelParams.fixta        =  *fixta;
			modelParams.nAsig0       =  *nAsig0;
			modelParams.nta0         =  *nta0;
			modelParams.Nsur         =  *Nsur;
			modelParams.LLinversion  =  *LLinversion;
			modelParams.forecast     =  *forecast;
			modelParams.r0  		 =  *r0;
			modelParams.Asig0        =  *Asig0;
			modelParams.ta0  		 =  *ta0;
			modelParams.Asig_min     =  *Asig_min;
			modelParams.Asig_max     =  *Asig_max;
			modelParams.ta_min       =  *ta_min;
			modelParams.ta_max       =  *ta_max;
			modelParams.tw   		 =  *tw;
			modelParams.fore_dt      =  *fore_dt;
			modelParams.Mc   		 =  *Mc;
			modelParams.Mag_main     =  *Mag_main;
			modelParams.Mc_source    =  *Mc_source;
			modelParams.dCFS         =  *dCFS;
			modelParams.DCFS_cap     =  *DCFS_cap;
			modelParams.dt 			 =  *dt;
			modelParams.dM  		 =  *dM;
			modelParams.xytoll       =  *xytoll;
			modelParams.ztoll        =  *ztoll;
			modelParams.border       =  *border;
			modelParams.res 		 =  *res;
			modelParams.gridresxy    =  *gridresxy;
			modelParams.gridresz     =  *gridresz;
			modelParams.smoothing    =  *smoothing;

		}

		//Map the BCast_Model_Parameters struct to MPI struct type
		MPI_Datatype CRS_MPI_BCast_Model_Parameters;
		int blocks_ModelParameters				[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Datatype types_ModelParameters		[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Aint addresses_ModelParameters		[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Aint displacements_ModelParameters	[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Aint baseAddress_ModelParameters;

		// Set blocks
		for(int i = 0; i < SIZE_BCAST_MODEL_PARAMETERS; ++i) {
			blocks_ModelParameters[i] = 1;
		}

		// Set types
		for(int i = 0; i < 8; ++i) {
			types_ModelParameters[i] = MPI_INT;
		}
		for(int i = 8; i < SIZE_BCAST_MODEL_PARAMETERS; ++i) {
			types_ModelParameters[i] = MPI_DOUBLE;
		}

		// Set addresses
		MPI_Address(&modelParams, &baseAddress_ModelParameters);

		MPI_Address(&(modelParams.fixr),        &addresses_ModelParameters[0]);
		MPI_Address(&(modelParams.fixAsig),     &addresses_ModelParameters[1]);
		MPI_Address(&(modelParams.fixta),       &addresses_ModelParameters[2]);
		MPI_Address(&(modelParams.nAsig0),      &addresses_ModelParameters[3]);
		MPI_Address(&(modelParams.nta0),        &addresses_ModelParameters[4]);
		MPI_Address(&(modelParams.Nsur),        &addresses_ModelParameters[5]);
		MPI_Address(&(modelParams.LLinversion), &addresses_ModelParameters[6]);
		MPI_Address(&(modelParams.forecast),    &addresses_ModelParameters[7]);

		MPI_Address(&(modelParams.r0), 			&addresses_ModelParameters[8]);
		MPI_Address(&(modelParams.Asig0),       &addresses_ModelParameters[9]);
		MPI_Address(&(modelParams.ta0), 		&addresses_ModelParameters[10]);
		MPI_Address(&(modelParams.Asig_min),    &addresses_ModelParameters[11]);
		MPI_Address(&(modelParams.Asig_max),    &addresses_ModelParameters[12]);
		MPI_Address(&(modelParams.ta_min),      &addresses_ModelParameters[13]);
		MPI_Address(&(modelParams.ta_max),      &addresses_ModelParameters[14]);
		MPI_Address(&(modelParams.tw),  		&addresses_ModelParameters[15]);
		MPI_Address(&(modelParams.fore_dt),     &addresses_ModelParameters[16]);
		MPI_Address(&(modelParams.Mc),  		&addresses_ModelParameters[17]);
		MPI_Address(&(modelParams.Mag_main),    &addresses_ModelParameters[18]);
		MPI_Address(&(modelParams.Mc_source),   &addresses_ModelParameters[19]);
		MPI_Address(&(modelParams.dCFS),        &addresses_ModelParameters[20]);
		MPI_Address(&(modelParams.DCFS_cap),    &addresses_ModelParameters[21]);
		MPI_Address(&(modelParams.dt),  		&addresses_ModelParameters[22]);
		MPI_Address(&(modelParams.dM),  		&addresses_ModelParameters[23]);
		MPI_Address(&(modelParams.xytoll),      &addresses_ModelParameters[24]);
		MPI_Address(&(modelParams.ztoll),       &addresses_ModelParameters[25]);
		MPI_Address(&(modelParams.border),      &addresses_ModelParameters[26]);
		MPI_Address(&(modelParams.res), 		&addresses_ModelParameters[27]);
		MPI_Address(&(modelParams.gridresxy),   &addresses_ModelParameters[28]);
		MPI_Address(&(modelParams.gridresz),    &addresses_ModelParameters[29]);
		MPI_Address(&(modelParams.smoothing),   &addresses_ModelParameters[30]);



		// Set displacements
		for(int i = 0; i < SIZE_BCAST_MODEL_PARAMETERS; ++i) {
			displacements_ModelParameters[i] = addresses_ModelParameters[i] - baseAddress_ModelParameters;
		}

		MPI_Type_struct(SIZE_BCAST_MODEL_PARAMETERS, blocks_ModelParameters,
						displacements_ModelParameters, types_ModelParameters,
						&CRS_MPI_BCast_Model_Parameters);

		MPI_Type_commit(&CRS_MPI_BCast_Model_Parameters);

		MPI_Bcast(&modelParams, 1, CRS_MPI_BCast_Model_Parameters, 0, MPI_COMM_WORLD);

		// For processes other than root, populate the scalars from the struct
		// received from root.
		if(procId != 0) {
			*fixr   		= modelParams.fixr;
			*fixAsig    	= modelParams.fixAsig;
			*fixta  		= modelParams.fixta;
			*nAsig0     	= modelParams.nAsig0;
			*nta0   		= modelParams.nta0;
			*Nsur   		= modelParams.Nsur;
			*LLinversion    = modelParams.LLinversion;
			*forecast       = modelParams.forecast;
			*r0     		= modelParams.r0;
			*Asig0  		= modelParams.Asig0;
			*ta0    		= modelParams.ta0;
			*Asig_min       = modelParams.Asig_min;
			*Asig_max       = modelParams.Asig_max;
			*ta_min         = modelParams.ta_min;
			*ta_max         = modelParams.ta_max;
			*tw     		= modelParams.tw;
			*fore_dt        = modelParams.fore_dt;
			*Mc     		= modelParams.Mc;
			*Mag_main       = modelParams.Mag_main;
			*Mc_source      = modelParams.Mc_source;
			*dCFS   		= modelParams.dCFS;
			*DCFS_cap       = modelParams.DCFS_cap;
			*dt     		= modelParams.dt;
			*dM     		= modelParams.dM;
			*xytoll         = modelParams.xytoll;
			*ztoll  		= modelParams.ztoll;
			*border         = modelParams.border;
			*res    		= modelParams.res;
			*gridresxy      = modelParams.gridresxy;
			*gridresz       = modelParams.gridresz;
			*smoothing      = modelParams.smoothing;

		}

		// Map struct flags to MPI struct type
		MPI_Datatype CRS_MPI_BCast_Flags;
		int blocks_Flags[BCAST_FLAGS_SIZE];
		MPI_Datatype types_Flags[BCAST_FLAGS_SIZE];
		MPI_Aint baseAddress_Flags;
		MPI_Aint addresses_Flags[BCAST_FLAGS_SIZE];
		MPI_Aint displacements_Flags[BCAST_FLAGS_SIZE];

		// Set blocks and types
		for(int i = 0; i < BCAST_FLAGS_SIZE; ++i) {
			blocks_Flags[i] = 1;
			types_Flags[i] = MPI_INT;
		}

		// Set addresses
		MPI_Address(flags, &baseAddress_Flags);

		MPI_Address(&(flags->err_recfault),      		&addresses_Flags[0]);
		MPI_Address(&(flags->err_gridpoints),    		&addresses_Flags[1]);
		MPI_Address(&(flags->OOPs),      		 		&addresses_Flags[2]);
		MPI_Address(&(flags->afterslip),         		&addresses_Flags[3]);
		MPI_Address(&(flags->splines),   		 		&addresses_Flags[4]);
		MPI_Address(&(flags->sources_all_iso),   		&addresses_Flags[5]);
		MPI_Address(&(flags->sources_without_focmec),   &addresses_Flags[6]);
		MPI_Address(&(flags->sample_all),        		&addresses_Flags[7]);

		// Set displacements
		for(int i = 0; i < BCAST_FLAGS_SIZE; ++i) {
			displacements_Flags[i] = addresses_Flags[i] - baseAddress_Flags;
		}

		MPI_Type_struct(BCAST_FLAGS_SIZE, blocks_Flags, displacements_Flags,
						types_Flags, &CRS_MPI_BCast_Flags);
		MPI_Type_commit(&CRS_MPI_BCast_Flags);

		MPI_Bcast(flags, 1, CRS_MPI_BCast_Flags, 0, MPI_COMM_WORLD);

		// Broadcast crust structure and related variables:
//		MPI_Bcast(&((*crst).lambda),	1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(&((*crst).mu), 		1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(&((*crst).fric), 		1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(&((*crst).skepton), 	1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(&((*crst).str0), 		1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(&((*crst).dip0), 		1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(&((*crst).rake0), 	1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(s, 	3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(st, 	3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(di, 	3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Bcast(regstress_mode, 		120, MPI_CHAR,   0, MPI_COMM_WORLD);

		MPI_Bcast(s,  3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(st, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(di, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Bcast(&(crst->lambda),	 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->mu), 		 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->fric), 	 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->skepton),  1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->str0[0]),  1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->dip0[0]),  1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->rake0[0]), 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(regstress_mode, 	 120, MPI_CHAR,   0, MPI_COMM_WORLD);

//		// FIXME: Fahad - For debugging purposes only ...
//		printf("\nProcId: %d -- (*crst).lambda: %f",  procId, (*crst).lambda);
//		printf("\nProcId: %d -- (*crst).mu: %f", 	  procId, (*crst).mu);
//		printf("\nProcId: %d -- (*crst).fric: %f", 	  procId, (*crst).fric);
//		printf("\nProcId: %d -- (*crst).skepton: %f", procId, (*crst).skepton);
//		printf("\nProcId: %d -- (*crst).str0[0]: %f", procId, (*crst).str0[0]);
//		printf("\nProcId: %d -- (*crst).dip0[0]: %f", procId, (*crst).dip0[0]);
//		printf("\nProcId: %d -- (*crst).rake0[0]: %f", procId, (*crst).rake0[0]);
//		for(int i = 0; i < 3; ++i) {
//			printf("\nProcId: %d -- s[%d]: %f",  procId, i, s[0]);
//			printf("\nProcId: %d -- st[%d]: %f", procId, i, st[0]);
//			printf("\nProcId: %d -- di[%d]: %f", procId, i, di[0]);
//		}
//		printf("\nProcId: %d -- regstress_mode: %s", procId, regstress_mode);
//
//		MPI_Barrier(MPI_COMM_WORLD);
//		error_quit("\n\nread_param.c -- Exiting at line 465 \n");

	#endif

	//calculate regional stress field:
	if (!(strcmp(regstress_mode,"oops"))) {
		prestress(1e6*s[0], 1e6*s[1], 1e6*s[2], (*crst).str0[0], (*crst).dip0[0], (*crst).rake0[0], 0.0,(*crst).fric, &((*crst).S));
	}
	else{
		if (!(strcmp(regstress_mode,"paxis"))) {
			for (int i=0; i<3; i++) s[i]*=1e6;
			(*crst).S=prestress_eigen(s, st, di);
		}
	}

	return 0;
}
