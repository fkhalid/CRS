/*
 * read_inputfile.c
 *
 *  Created on: Nov 21, 2013
 *      Author: camcat
 */

#include "read_inputfile.h"

int read_inputfile(char *input_fname, char *outname, char *reftime_str, char *crust_file, char *fore_template,
		char *catname, char ***focmeccat, char *background_rate_file, char *slipmodelfile, char *afterslipmodelfile,
		char *model_parameters_file, char *Logfile, int *extraoutput, struct tm *reftime,
		double *Tstart, double *Tend, long *seed, char *cmb_format, int *num_fm){

	/* input: file name input_fname
	 *
	 * output:
	 * 		outname: output file
	 * 		reftime_str: string containing the reference time (issue time)
	 * 		reftime: structure containing the reference time (issue time)
	 * 		Tstart, Tend: time in days from reftime
	 * 		crust_file: input file for farfalle code (general file with info about the crust etc).
	 * 		fore_template: forecast template
	 * 		catname: catalog
	 * 		focmeccat: catalog of focal mechanisms
	 * 		background_rate_file: file containing background seismicity model
	 * 		slipmodefile, afterslipmodelfile: files containing a list of slip models/afterslipmodel snapshots
	 *
	 * NB: all pointers will be ignored if NULL is passed; otherwise, char* should already be initialized.
	 */

	setenv("TZ", "UTC", 1);

	FILE *fin;
	int Nchar=1000;
	char line[Nchar], listfocmeccat[Nchar];
	char *key, *value;
	int NP=17, i, err=0;
	struct tm times;
	int value_found[NP], listfm=0, nofm=0;
	char comment[]="#", comm=comment[0];

	for (int n=0; n<NP; n++) value_found[n]=0;

	char *keys[]={
	/*0*/	"IssueDate",
	/*1*/	"ForecastStartDate", \
	/*2*/	"ForecastEndDate",	\
	/*3*/	"OutputForecastFile", \
	/*4*/	"InputCatalogFile",	\
	/*5*/	"InputCatalogFocMecFile",	\
	/*6*/	"InputListCatalogFocMecFile",	\
	/*7*/	"ForecastTemplate", \
	/*8*/	"InputCoulombFile", \
	/*9*/	"InputListSlipModels", \
	/*10*/	"InputListAfterslipModels", \
	/*11*/	"InputBackgroundRateFile",\
	/*12*/	"InputModelParametersFile",\
	/*13*/	"RandomSeedValue",\
	/*14*/	"Logfile",\
	/*15*/	"ExtendedOutput",\
	/*16*/	"CmbFormat"
	};


	if((fin = fopen(input_fname, "r"))==NULL) {
		if (verbose_level>1) fprintf(stderr, "Error read_input: unable to open input file %s.\n", input_fname);
		return 1;
	}

	while(!feof(fin)) {
		fgets(line,Nchar,fin);
		if (line[0]==comm) continue;
		if (ferror(fin)) {
			if (verbose_level>1) error_quit("Error reading input data using fgets!\n");
			return 1;
		}
		key=strtok(line,"=");
		value=strtok(NULL,"=");
		if (!value) continue;
		i=0;
		while (i<NP && strcmp(key,keys[i])) i++;
		if (i>=NP){
			if (verbose_level>0) fprintf(stderr, "Error read_input: parameter \" %s\" in file \"%s\" not recognized.\n", key, input_fname);
			continue;
		}

		value_found[i]=1;

		switch(i){
			case 0:
				sscanf(value,"%s",reftime_str);
				sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
				times.tm_year-=1900;
				times.tm_mon-=1;
				times.tm_isdst=0;
				if (reftime) *reftime=times;
				break;
			case 1:
				sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
				times.tm_year-=1900;
				times.tm_mon-=1;
				times.tm_isdst=0;
				if (Tstart) *Tstart=difftime(mktime(&times),mktime(reftime))*SEC2DAY;
				break;
			case 2:
				sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
				times.tm_year-=1900;
				times.tm_mon-=1;
				times.tm_isdst=0;
				if (Tend) *Tend=difftime(mktime(&times),mktime(reftime))*SEC2DAY;
				break;
			case 3:
				if (outname) sscanf(value,"%s",outname);
				break;
			case 4:
				if (catname) sscanf(value,"%s",catname);
				break;
			case 5:
				if (focmeccat) {
					*focmeccat= malloc(sizeof(char*));
					(*focmeccat)[0]=malloc(120*sizeof(char));
					sscanf(value,"%s",(*focmeccat)[0]);
				}
				if (num_fm) *num_fm=1;
				break;
			case 6:
				if (focmeccat) sscanf(value,"%s",listfocmeccat);
				listfm=1;
				break;
			case 7:
				if (fore_template) sscanf(value,"%s",fore_template);
				break;
			case 8:
				if (crust_file) sscanf(value,"%s",crust_file);
				break;
			case 9:
				if (slipmodelfile) sscanf(value,"%s",slipmodelfile);
				break;
			case 10:
				if (afterslipmodelfile) sscanf(value,"%s",afterslipmodelfile);
				break;
			case 11:
				if (background_rate_file) sscanf(value,"%s",background_rate_file);
				break;
			case 12:
				if (model_parameters_file) sscanf(value,"%s",model_parameters_file);
				break;
			case 13:
				if (seed) sscanf(value,"%ld",seed);
				break;
			case 14:
				if (Logfile) sscanf(value,"%s",Logfile);
				break;
			case 15:
				if (extraoutput) sscanf(value,"%d",extraoutput);
				break;
			case 16:
				if (cmb_format) sscanf(value,"%s",cmb_format);
				break;
		}
    }

	fclose(fin);

	if (focmeccat && listfm){
		err=read_slipformecfiles(listfocmeccat, focmeccat, num_fm);
		if (err) {
			if (verbose_level) printf("Error: could not read file %s.\n", listfocmeccat);
			if (flog) fprintf(flog,"Error: could not read file %s.\n", listfocmeccat);
		}
	}

	nofm=0;
	for (int n=0; n<NP; n++) {
		if (!value_found[n]) {
			switch (n){
			case 3:
				if (verbose_level>2) printf("Warning: parameter %s not given in %s -> will use output/forecast.\n", keys[n], input_fname);
				if (outname) strcpy(outname,"output/forecast");
				break;
			case 5:
				if (nofm){
					if (verbose_level>2) printf("Warning: parameters %s, %s not given in %s.\n", keys[5], keys[6], input_fname);
					if (focmeccat) *focmeccat=NULL;
				}
				else nofm=1;
				break;
			case 6:
				if (nofm){
					if (verbose_level>2) printf("Warning: parameters %s, %s not given in %s.\n", keys[5], keys[6], input_fname);
					if (focmeccat) *focmeccat=NULL;
				}
				else nofm=1;
				break;
			case 9:
				if (verbose_level>2) printf("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
				if (slipmodelfile) strcpy(slipmodelfile,"");
				break;
			case 10:
				if (verbose_level>2) printf("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
				if (afterslipmodelfile) strcpy(afterslipmodelfile,"");
				break;
			case 11:
				if (verbose_level>2) printf("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
				if (background_rate_file) strcpy(background_rate_file,"");
				break;
			case 12:
				if (Logfile) strcpy(Logfile,"");
				break;
			case 13:
				if (extraoutput) extraoutput=0;
				break;
			default:
				if (verbose_level) printf("Error: parameter %s not given in %s.\n", keys[n], input_fname);
				if (flog) fprintf(flog, "Error: parameter %s not given in %s.\n", keys[n], input_fname);
				return 1;
				break;
			}
		}
	}

	return 0;
}

int read_slipformecfiles(char *inputfile, char ***listfiles, int *nfiles){

	int Nchar=1000;
	char line[Nchar];
	char comment[]="#", comm=comment[0];
	FILE *fin;

	if (!(fin=fopen(inputfile,"r"))) {
		if (verbose_level) printf("Error: can not open file %s (read_slipformecfiles), Exiting.\n", inputfile);
		if (flog) fprintf(flog, "Error: can not open file %s (read_slipformecfiles), Exiting.\n", inputfile);
		*listfiles=NULL;
		return 1;
	}
	else {
		line[0]=comm;
		while (line[0]==comm)fgets(line,Nchar,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data (file: %s) using fgets!\n", inputfile);
		sscanf(line,"%d", nfiles);

		*listfiles = malloc((*nfiles)*sizeof(char*));
		for (int nn=0; nn<(*nfiles); nn++) {
			(*listfiles)[nn] = malloc(120 * sizeof(char));
			line[0]=comm;
			while (line[0]==comm)fgets(line,Nchar,fin);
			if (ferror(fin)) {
				fprintf(stderr, "ERROR reading input data (file: %s) using fgets!\n", inputfile);
				return 1;
			}
			sscanf(line,"%s", (*listfiles)[nn]);
		}
	}
	return 0;
}

int read_listslipmodel(char *input_fname, struct tm reftime, struct slipmodels_list *allslipmodels, double res, int is_afterslip){

	FILE *fin;
	int Nchar=1000;
	char line[Nchar];
	char time_str[50];
	char comment[]="#", comm=comment[0];
	double t0;
	struct tm times;
	int Nm0, nsm, no_slipmod;

	if((fin = fopen(input_fname, "r"))==NULL) {
		if (verbose_level>1) fprintf(stderr, "Warning read_input: no slip model file found (read_listslipmodel).\n");
		if (flog) fprintf(flog, "\nWarning read_input: no slip model file found (read_listslipmodel).\n");
		(*allslipmodels).NSM=0;
		(*allslipmodels).is_afterslip=is_afterslip;
		(*allslipmodels).tmain= NULL;
		(*allslipmodels).mmain= NULL;
		(*allslipmodels).disc= NULL;
		(*allslipmodels).Nfaults=NULL;
		(*allslipmodels).no_slipmodels=NULL;
		return 1;
	}

		else {
		line[0]=comm;
		while (line[0]==comm)fgets(line,Nchar,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data (file: %s) using fgets!\n", input_fname);
		sscanf(line,"%d", &Nm0);
		fgets(line,Nchar,fin);
		sscanf(line,"%s %d", cmb_format, &((*allslipmodels).constant_geometry));
		if (is_afterslip) {
			fgets(line,Nchar,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
			sscanf(line, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
			times.tm_year-=1900;
			times.tm_mon-=1;
			times.tm_isdst=0;
			t0=difftime(mktime(&times),mktime(&reftime))*SEC2DAY;
		}
		(*allslipmodels).NSM=Nm0;
		(*allslipmodels).is_afterslip=is_afterslip;
		(*allslipmodels).tmain=dvector(0,Nm0-1);
		(*allslipmodels).mmain= (is_afterslip)? NULL : dvector(0,Nm0-1);
		if (is_afterslip){
			(*allslipmodels).disc=dvector(0,0);
			(*allslipmodels).disc[0]=res;
			(*allslipmodels).Nfaults=ivector(0,0);
			(*allslipmodels).no_slipmodels=ivector(0,0);
			(*allslipmodels).no_slipmodels[0]=1;
		}
		else {
			(*allslipmodels).disc=dvector(0,Nm0-1);
			(*allslipmodels).Nfaults=ivector(0,Nm0-1);
			(*allslipmodels).no_slipmodels=ivector(0,Nm0-1);
		}
		(*allslipmodels).slipmodels = malloc(Nm0*sizeof(char*));
		nsm=0;
		for (int nn=0; nn<Nm0; nn++) {
			if (is_afterslip){
				(*allslipmodels).slipmodels[nn] = malloc(120 * sizeof(char));
				(*allslipmodels).Nfaults[nn]=1;	//actual value found later.
				fgets(line,Nchar,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%lf %s", (*allslipmodels).tmain+nn, (*allslipmodels).slipmodels[nn]);
				(*allslipmodels).tmain[nn]+=t0;
			}
			else{
				(*allslipmodels).disc[nn] = res;
				fgets(line,Nchar,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%s %lf %d", time_str, (*allslipmodels).mmain+nn, &no_slipmod);
				sscanf(time_str, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
				times.tm_year-=1900;
				times.tm_mon-=1;
				times.tm_isdst=0;
				 (*allslipmodels).tmain[nn]=difftime(mktime(&times),mktime(&reftime))*SEC2DAY;

				 (*allslipmodels).no_slipmodels[nn]=no_slipmod;
				 if (nsm+1+no_slipmod>Nm0) (*allslipmodels).slipmodels=realloc((*allslipmodels).slipmodels, (nsm+1+no_slipmod) * sizeof(char*));

				 for (int n=1; n<=no_slipmod; n++){
					(*allslipmodels).slipmodels[nsm] = malloc(120 * sizeof(char));
					fgets(line,Nchar,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
					sscanf(line,"%s", (*allslipmodels).slipmodels[nsm]);
					nsm++;

				}
			}
		}
		fclose(fin);
	}

	if (flog) {
		nsm=0;
		if (is_afterslip) fprintf(flog, "\nAfterslip input file: %s.\n", input_fname);
		else fprintf(flog, "\nSlip input file: %s.\n", input_fname);
		fprintf(flog, "%d %s slip models:\n", (*allslipmodels).NSM, is_afterslip? "after" : "");
		for (int m=0; m<(*allslipmodels).NSM; m++){
			if (is_afterslip){
				if (m==0) fprintf(flog, "\t time \t name\n");
				fprintf(flog, "\t%.2lf\t%s\n", (*allslipmodels).tmain[m], (*allslipmodels).slipmodels[m]);
			}
			else{
				if (m==0) fprintf(flog, "\t time \t mag \t name\n");
				for (int n=1; n<=(*allslipmodels).no_slipmodels[m]; n++){
					fprintf(flog, "\t%.2lf\t%.2lf\t%s\n", (*allslipmodels).tmain[m], (*allslipmodels).mmain[m], (*allslipmodels).slipmodels[nsm]);
					nsm++;
				}
			}
		}
		fflush(flog);
	}

	return 0;
}

