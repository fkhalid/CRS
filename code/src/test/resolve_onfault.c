/*
 * resolve_onfault.c
 *
 *  Created on: Jun 24, 2014
 *      Author: camcat
 *
 *  Resolve DCFS from a slip model on a given location and rec. fault orientation.
 *
 */


#include <stdio.h>
#include <string.h>

#include "../defines.h"
//#include "../inp_out/read_crust.h"
//#include "../inp_out/read_eqkfm.h"
//#include "../okada/okadaDCFS.h"
//#include "../util/error.h"
//#include "../util/moreutil.h"
//#include "../util/nrutil.h"

int resolve_onfault(char **argv){

	int err=0;
	int nchar=200;
	char msg[nchar], infile[nchar], slipmodelfile[nchar];
	double lat, lon, dep, strike, dip, rake;
	int optrake, NF;
	struct pscmp DCFS;
	struct eqkfm *eqkfm;
	struct crust crst;


	printf("Input file: %s\n", argv[1]);
	fflush(stdout);
	sscanf(argv[1],"%s", infile);

	sscanf("testinput.txt","%s", infile);

	err=read_inputfile2(infile, slipmodelfile, &lat, &lon, &dep, &strike, &dip, &rake, &optrake);
	if (err) {
		sprintf(msg,"Error reading input file %s.\n", infile);
		error_quit(msg);
	}
	init_crst(&crst);
	err=read_pscmp_crust(slipmodelfile, &crst);
	crst.nofmzones=1;
	crst.x=dvector(1,1);
	crst.y=dvector(1,1);
	crst.depth=dvector(1,1);
	crst.x[1]=0.0;
	crst.y[1]=0.0;
	crst.depth[1]=dep;

	err+=read_pscmp_eqkfm(slipmodelfile, &eqkfm, &NF);

	//--------------------------------------------------------------------//
//	struct crust crst2;
//	double lat0=38, lon0=141.5, dep0=10.0;
//	int Nlat=50, Nlon=50, Ndep=1, NP=Nlat*Nlon*Ndep;
//	double *latgrid, *longrid, *depgrid;
//	int *allp;
//	double dlat=8.0, dlon=7.0, ddep=10.0;
//
//	latgrid=dvector(1,NP);
//	longrid=dvector(1,NP);
//	depgrid=dvector(1,NP);
//	allp=ivector(1,NP);
//
//	for (int i=1; i<=NP; i++) allp[i]=i;
//
//	for (int i=1; i<=Nlat; i++){
//		for (int j=1; j<=Nlon; j++){
//			for (int k=1; k<=Ndep; k++){
//				latgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lat0-0.5*dlat+(dlat/Nlat)*(i-0.5);
//				longrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lon0-0.5*dlon+(dlon/Nlon)*(j-0.5);
//				depgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= (ddep/Ndep)*(k-0.5);
//			}
//		}
//	}
//
//	init_crst(&crst2);
//	read_pscmp_crust(slipmodelfile,&crst2);
//
//	crst2.N_allP=NP;
//	crst2.lat=latgrid;
//	crst2.lon=longrid;
//	crst2.depth=depgrid;
//	crst2.x=dvector(1,NP);
//	crst2.y=dvector(1,NP);
//
//	DCFS.nsel=NP;
//	DCFS.which_pts=ivector(1,NP);
//	DCFS.cmb=dvector(1,NP);
//	DCFS.S=d3tensor(1,NP,1,3,1,3);
//	for (int i=1; i<=NP; i++) {
//		DCFS.which_pts[i]=i;
//		latlon2localcartesian(crst2.lat[i], crst2.lon[i], lat, lon, &(crst2.y[i]), &(crst2.x[i]));
//
//	}
	//--------------------------------------------------------------------//


	if (err) {
		sprintf(msg,"Error reading input file %s.\n", slipmodelfile);
		error_quit(msg);
	}
	for (int i=0; i<NF; i++ ){
		eqkfm[i].nsel=1;
		latlon2localcartesian(eqkfm[i].lat, eqkfm[i].lon, lat, lon, &(eqkfm[i].y), &(eqkfm[i].x));
	}

	DCFS.nsel=1;
	DCFS.which_pts=ivector(1,1);
	DCFS.cmb=dvector(1,1);
	DCFS.S=d3tensor(1,1,1,3,1,3);
	DCFS.which_pts[1]=1;

	okadaDCFS(DCFS, eqkfm, NF, crst, NULL, NULL, 1);
	resolve_DCFS(DCFS, crst, &strike, &dip, &rake, optrake);

	printf("\n -------------------------------------\n           DCFS=%.3e MPa. \n -------------------------------------\n ", DCFS.cmb[1]*1e-6);

//	print_grid("dcfs2.dat", DCFS, crst2, DCFS.cmb);


	return err;
}


int read_inputfile2(char *input_fname, char *slipmodelfile, double *lat, double *lon, double *dep, double *str, double *dip, double *rake, int *optrake){
	/* input: file name input_fname
	 *
	 * output:
	 * 		slipmodefile: file containing slip model, in pscmp format.
	 *		lon, lat, dep, str, dip, rake: location and orientation of receiver point.
	 */


	FILE *fin;
	int Nchar=1000;
	char line[Nchar];
	char *key, *value;
	int NP=7, i, err=0;
	int value_found[NP];
	char comment[]="#", comm=comment[0];

	for (int n=0; n<NP; n++) value_found[n]=0;

	char *keys[]={
	/*0*/	"slip_model",
	/*1*/	"lat", \
	/*2*/	"lon",	\
	/*3*/	"dep", \
	/*4*/	"strike",	\
	/*5*/	"dip",	\
	/*6*/	"rake"	//optional (optrake=1 if not set).
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
				sscanf(value,"%s",slipmodelfile);
				break;
			case 1:
				sscanf(value, "%lf", lat);
				break;
			case 2:
				sscanf(value, "%lf", lon);
				break;
			case 3:
				sscanf(value, "%lf", dep);
				break;
			case 4:
				sscanf(value, "%lf", str);
				break;
			case 5:
				sscanf(value, "%lf", dip);
				break;
			case 6:
				sscanf(value, "%lf", rake);
				*optrake=0;
				break;
		}
    }

	fclose(fin);

	for (int n=0; n<NP; n++) {
		if (!value_found[n]) {
			switch (n){
			case 6:
				printf("Warning: rake not given, will calculate optimal rake in total stress field (background + DCFS).\n");
				*optrake=1;
				break;
			default:
				printf("Error: parameter %s not given in %s.\n", keys[n], input_fname);
				return 1;
				break;
			}
		}
	}

	return err;
}
