/*
 * read_param.c
 *
 *  Created on: Dec 20, 2013
 *      Author: camcat
 */

#include <stdio.h>
#include <time.h>

#include "../defines.h"

int read_modelparmeters(char * modelparametersfile, struct tm reftime, int *N_min_events, int * fixr, int *fixAsig, int *fixta, double *r0, double *Asig0, double *ta0,
		double *Asig_min, double *Asig_max, double *ta_min, double *ta_max, int *nAsig0, int *nta0, double *tstartLL, double *extra_time, double *tw, double *fore_dt, double *t_back,
		int *Nsur, int *Nslipmod, struct flags *flags, double *Hurst, double *Mc_source, int *use_bg_rate, double *Mc, double *Mag_main, double *DCFS_cap, int *gridPMax,
		double *dt, double *dM, double *xytoll, double *ztoll, double *border, double *res, double *gridresxy, double *gridresz, double *smoothing, int *LLinversion, int *forecast){

	FILE * fin;
	char comment[]="#", comm=comment[0];
	int Nchar_long=500;
	char line[Nchar_long];
	struct tm times;
	int aftershock_mode;

	if (!(fin=fopen(modelparametersfile,"r"))){
		if (verbose_level) printf("Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);
		if (flog) fprintf(flog, "Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);
		return 1;
	}
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
	sscanf(line,"%d %d %lf", &((*flags).err_slipmodel), &((*flags).err_afterslipmodel), Hurst);
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

	return 0;
}
