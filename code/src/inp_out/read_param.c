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

int read_modelparameters(char * modelparametersfile, struct tm reftime, int *N_min_events,
						int *fixr, int *fixAsig, int *fixta, double *r0, double *Asig0,
						double *ta0, double *Asig_min, double *Asig_max, double *ta_min,
						double *ta_max, int *nAsig0, int *nta0, double *tstartLL,
						double *extra_time, double *tw, double *fore_dt, double *t_back,
						int *Nsur, int *Nslipmod, struct flags *flags,
						double *Mc_source, int *use_bg_rate, double *Mc, double *Mag_main,
						double *DCFS_cap, int *gridPMax, double *dt, double *dM,
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
	struct tm times;
	int aftershock_mode;

	// [Fahad] If there is a file error, only root will know about it.
	//		   So it is important that the error is broadcast to all
	//		   processes, and they all return from the function with
	//		   the same error code.
	if(procId == 0) {
		fin = fopen(modelparametersfile,"r");
		if(fin == NULL) {
			if (verbose_level) printf("Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);
			if (flog) fprintf(flog, "Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);

			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	if(procId == 0) {
		sprintf(comment,"#");
		comm=comment[0];
		line[0]=comm;
		while (line[0]==comm)fgets(line,Nchar_long,fin);
		sscanf(line,"%d", N_min_events);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf", fixr, r0);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf %lf %lf %d", fixAsig, Asig0, Asig_min, Asig_max, nAsig0);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf %lf %lf %d", fixta, ta0, ta_min, ta_max, nta0);
		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
		times.tm_year-=1900;
		times.tm_mon-=1;
		times.tm_isdst=0;
		*tstartLL=difftime(mktime(&times),mktime(&reftime))*SEC2DAY;
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", tw);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", extra_time);
		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%lf", fore_dt);
		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%d %d", Nsur, Nslipmod);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		//sscanf(line,"%d %d %lf", &((*flags).err_slipmodel), &((*flags).err_afterslipmodel), Hurst);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", &((*flags).err_recfault));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", &((*flags).err_gridpoints));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", &((*flags).OOPs));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", &((*flags).afterslip));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf %d", &((*flags).aftershocks), Mc_source, &aftershock_mode);
		switch (aftershock_mode){
			case 0:
				(*flags).only_aftershocks_withfm=0;
				(*flags).full_field=0;
				(*flags).aftershocks_fixedmec=0;
				break;
			case 1:
				(*flags).only_aftershocks_withfm=0;
				(*flags).full_field=1;
				(*flags).aftershocks_fixedmec=0;
				break;
			case 2:
				(*flags).only_aftershocks_withfm=0;
				(*flags).full_field=2;
				(*flags).aftershocks_fixedmec=1;
				break;
			case 3:
				(*flags).only_aftershocks_withfm=0;
				(*flags).full_field=2;
				(*flags).aftershocks_fixedmec=0;
				break;
			case 4:
				(*flags).only_aftershocks_withfm=1;
				(*flags).full_field=2;
				(*flags).aftershocks_fixedmec=0;
				break;
			default:
				break;
		}

		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", use_bg_rate);
		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%lf", Mc);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", Mag_main);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", DCFS_cap);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", gridPMax);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", dt);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", dM);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", xytoll);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", ztoll);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", border);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", res);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", gridresxy);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", gridresz);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", smoothing);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
		times.tm_year-=1900;
		times.tm_mon-=1;
		times.tm_isdst=0;
		*t_back=difftime(mktime(&times),mktime(&reftime))*SEC2DAY;
		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%d", LLinversion);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", forecast);

		fclose(fin);
	}

	// TODO: [Fahad] Look into how much of the following code can be moved to
	//		 a separate function ...
	#ifdef _CRS_MPI
		// Copy scalars to the BCast_Model_Parameters struct
		struct BCast_Model_Parameters modelParams;

		if(procId == 0) {
			modelParams.N_min_events = *N_min_events;
			modelParams.fixr 		 = *fixr;
			modelParams.fixAsig 	 = *fixAsig;
			modelParams.fixta 		 = *fixta;
			modelParams.nAsig0 		 = *nAsig0;
			modelParams.nta0 		 = *nta0;
			modelParams.Nsur 		 = *Nsur;
			modelParams.Nslipmod 	 = *Nslipmod;
			modelParams.use_bg_rate  = *use_bg_rate;
			modelParams.gridPMax 	 = *gridPMax;
			modelParams.LLinversion  = *LLinversion;
			modelParams.forecast 	 = *forecast;
			modelParams.r0 			 = *r0;
			modelParams.Asig0 		 = *Asig0;
			modelParams.ta0 		 = *ta0;
			modelParams.Asig_min 	 = *Asig_min;
			modelParams.Asig_max 	 = *Asig_max;
			modelParams.ta_min 		 = *ta_min;
			modelParams.ta_max 		 = *ta_max;
			modelParams.tstartLL 	 = *tstartLL;
			modelParams.extra_time 	 = *extra_time;
			modelParams.tw 			 = *tw;
			modelParams.fore_dt 	 = *fore_dt;
			modelParams.t_back 		 = *t_back;
			modelParams.Hurst 		 = *Hurst;
			modelParams.Mc_source 	 = *Mc_source;
			modelParams.Mc 			 = *Mc;
			modelParams.Mag_main 	 = *Mag_main;
			modelParams.DCFS_cap 	 = *DCFS_cap;
			modelParams.dt 			 = *dt;
			modelParams.dM 			 = *dM;
			modelParams.xytoll 		 = *xytoll;
			modelParams.ztoll 		 = *ztoll;
			modelParams.border 		 = *border;
			modelParams.res 		 = *res;
			modelParams.gridresxy 	 = *gridresxy;
			modelParams.gridresz 	 = *gridresz;
			modelParams.smoothing 	 = *smoothing;
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
		for(int i = 0; i < 12; ++i) {
			types_ModelParameters[i] = MPI_INT;
		}
		for(int i = 12; i < SIZE_BCAST_MODEL_PARAMETERS; ++i) {
			types_ModelParameters[i] = MPI_DOUBLE;
		}

		// Set addresses
		MPI_Address(&modelParams, &baseAddress_ModelParameters);

		MPI_Address(&(modelParams.N_min_events),	&addresses_ModelParameters[0]);
		MPI_Address(&(modelParams.fixr), 			&addresses_ModelParameters[1]);
		MPI_Address(&(modelParams.fixAsig), 		&addresses_ModelParameters[2]);
		MPI_Address(&(modelParams.fixta), 			&addresses_ModelParameters[3]);
		MPI_Address(&(modelParams.nAsig0), 			&addresses_ModelParameters[4]);
		MPI_Address(&(modelParams.nta0), 			&addresses_ModelParameters[5]);
		MPI_Address(&(modelParams.Nsur), 			&addresses_ModelParameters[6]);
		MPI_Address(&(modelParams.Nslipmod), 		&addresses_ModelParameters[7]);
		MPI_Address(&(modelParams.use_bg_rate), 	&addresses_ModelParameters[8]);
		MPI_Address(&(modelParams.gridPMax), 		&addresses_ModelParameters[9]);
		MPI_Address(&(modelParams.LLinversion), 	&addresses_ModelParameters[10]);
		MPI_Address(&(modelParams.forecast), 		&addresses_ModelParameters[11]);
		MPI_Address(&(modelParams.r0), 				&addresses_ModelParameters[12]);
		MPI_Address(&(modelParams.Asig0), 			&addresses_ModelParameters[13]);
		MPI_Address(&(modelParams.ta0), 			&addresses_ModelParameters[14]);
		MPI_Address(&(modelParams.Asig_min), 		&addresses_ModelParameters[15]);
		MPI_Address(&(modelParams.Asig_max), 		&addresses_ModelParameters[16]);
		MPI_Address(&(modelParams.ta_min), 			&addresses_ModelParameters[17]);
		MPI_Address(&(modelParams.ta_max), 			&addresses_ModelParameters[18]);
		MPI_Address(&(modelParams.tstartLL), 		&addresses_ModelParameters[19]);
		MPI_Address(&(modelParams.extra_time), 		&addresses_ModelParameters[20]);
		MPI_Address(&(modelParams.tw), 				&addresses_ModelParameters[21]);
		MPI_Address(&(modelParams.fore_dt), 		&addresses_ModelParameters[22]);
		MPI_Address(&(modelParams.t_back), 			&addresses_ModelParameters[23]);
		MPI_Address(&(modelParams.Hurst), 			&addresses_ModelParameters[24]);
		MPI_Address(&(modelParams.Mc_source), 		&addresses_ModelParameters[25]);
		MPI_Address(&(modelParams.Mc), 				&addresses_ModelParameters[26]);
		MPI_Address(&(modelParams.Mag_main), 		&addresses_ModelParameters[27]);
		MPI_Address(&(modelParams.DCFS_cap), 		&addresses_ModelParameters[28]);
		MPI_Address(&(modelParams.dt), 				&addresses_ModelParameters[29]);
		MPI_Address(&(modelParams.dM), 				&addresses_ModelParameters[30]);
		MPI_Address(&(modelParams.xytoll),	 		&addresses_ModelParameters[31]);
		MPI_Address(&(modelParams.ztoll), 			&addresses_ModelParameters[32]);
		MPI_Address(&(modelParams.border), 			&addresses_ModelParameters[33]);
		MPI_Address(&(modelParams.res), 			&addresses_ModelParameters[34]);
		MPI_Address(&(modelParams.gridresxy), 		&addresses_ModelParameters[35]);
		MPI_Address(&(modelParams.gridresz), 		&addresses_ModelParameters[36]);
		MPI_Address(&(modelParams.smoothing), 		&addresses_ModelParameters[37]);

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
			*N_min_events 	= modelParams.N_min_events;
			*fixr 			= modelParams.fixr;
			*fixAsig 		= modelParams.fixAsig;
			*fixta 			= modelParams.fixta;
			*nAsig0 		= modelParams.nAsig0;
			*nta0 			= modelParams.nta0;
			*Nsur 			= modelParams.Nsur;
			*Nslipmod 		= modelParams.Nslipmod;
			*use_bg_rate 	= modelParams.use_bg_rate;
			*gridPMax 		= modelParams.gridPMax;
			*LLinversion 	= modelParams.LLinversion;
			*forecast 		= modelParams.forecast;
			*r0 			= modelParams.r0;
			*Asig0 			= modelParams.Asig0;
			*ta0 			= modelParams.ta0;
			*Asig_min 		= modelParams.Asig_min;
			*Asig_max 		= modelParams.Asig_max;
			*ta_min 		= modelParams.ta_min;
			*ta_max 		= modelParams.ta_max;
			*tstartLL 		= modelParams.tstartLL;
			*extra_time 	= modelParams.extra_time;
			*tw 			= modelParams.tw;
			*fore_dt 		= modelParams.fore_dt;
			*t_back 		= modelParams.t_back;
			*Hurst 			= modelParams.Hurst;
			*Mc_source 		= modelParams.Mc_source;
			*Mc 			= modelParams.Mc;
			*Mag_main 		= modelParams.Mag_main;
			*DCFS_cap 		= modelParams.DCFS_cap;
			*dt 			= modelParams.dt;
			*dM 			= modelParams.dM;
			*xytoll 		= modelParams.xytoll;
			*ztoll 			= modelParams.ztoll;
			*border 		= modelParams.border;
			*res 			= modelParams.res;
			*gridresxy 		= modelParams.gridresxy;
			*gridresz 		= modelParams.gridresz;
			*smoothing 		= modelParams.smoothing;
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

		MPI_Address(&(flags->err_recfault), 			&addresses_Flags[0]);
		MPI_Address(&(flags->err_slipmodel), 			&addresses_Flags[1]);
		MPI_Address(&(flags->err_afterslipmodel), 		&addresses_Flags[2]);
		MPI_Address(&(flags->err_gridpoints), 			&addresses_Flags[3]);
		MPI_Address(&(flags->OOPs), 					&addresses_Flags[4]);
		MPI_Address(&(flags->afterslip), 				&addresses_Flags[5]);
		MPI_Address(&(flags->splines), 					&addresses_Flags[6]);
		MPI_Address(&(flags->aftershocks), 				&addresses_Flags[7]);
		MPI_Address(&(flags->only_aftershocks_withfm), 	&addresses_Flags[8]);
		MPI_Address(&(flags->full_field), 				&addresses_Flags[9]);
		MPI_Address(&(flags->aftershocks_fixedmec), 	&addresses_Flags[10]);
		MPI_Address(&(flags->aftershocks_mode), 		&addresses_Flags[11]);
		MPI_Address(&(flags->new_slipmodel), 			&addresses_Flags[12]);
		MPI_Address(&(flags->sample_all), 				&addresses_Flags[13]);

		// Set displacements
		for(int i = 0; i < BCAST_FLAGS_SIZE; ++i) {
			displacements_Flags[i] = addresses_Flags[i] - baseAddress_Flags;
		}

		MPI_Type_struct(BCAST_FLAGS_SIZE, blocks_Flags, displacements_Flags,
						types_Flags, &CRS_MPI_BCast_Flags);
		MPI_Type_commit(&CRS_MPI_BCast_Flags);

		MPI_Bcast(flags, 1, CRS_MPI_BCast_Flags, 0, MPI_COMM_WORLD);
	#endif

	return 0;
}
