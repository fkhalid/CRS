/*
 * read_crust.c
 *
 *  Created on: Jan 16, 2013
 *      Author: camcat
 */

#include "read_crust.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

//TODO: if more slip models are contained, make sure they have consistent geometry (they must cover all sampling points).

int read_crust(char *fname, char *fnametemplate, char *focmecgridfile, struct crust *crst, double resxy, double resz){
/*
 * Read crust master file into crst structure.
 *
 * input:
 * 	fname= file containing info about the crust (pscmp or farfalle format).
 * 	fnametemplate= grid file (CSEP format).
 * 	focmecgridfile= grid file containing foc. planes for each grid point.
 *  resxy, resz= resired grid resolution (for calculations);
 *
 * output:
 * 	crst= structure containing info about the domain;
 *
 */

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int err=0, err1=0;
	int NG, ind, NLat, NLon, Nd;
	int no_subpointsx, no_subpointsy, no_subpointsz;
	int no_magbins;
	int is_refined=0;
	double mag1, mag2;
	double dx, dy, dAeq;
	double lat0, lat1,lon0, lon1, d0, d1;
	double *strtmp=0, *diptmp=0, *raketmp=0;

	if(procId == 0) {
		if (verbose_level>0) printf("Loading model setup...");
		if (flog) fprintf(flog,"\nEntering read_crust...\n");
	}

	double *olats, *olons, *odeps;

	//--------------read general crust information:-------------//

	init_crst(crst);
	if (!(strcmp(cmb_format,"farfalle"))) {
		err=read_farfalle_crust(fname, crst);
		if(procId == 0) {
			if (flog) fprintf(flog,"reading farfalle format (file %s)\n",fname);
		}
	}
	else {
		if (!(strcmp(cmb_format,"pscmp"))) {
			err=read_pscmp_crust(fname, crst);
			if(procId == 0) {
				if (flog) fprintf(flog,"reading pscmp format (file %s)\n",fname);
			}
		}
		else {
			if(procId == 0) {
				if (flog) fprintf(flog,"Unknown format: %s.\n", cmb_format);
			}
		}
	}

	if (err) {
		if(procId == 0) {
			if (flog) fprintf(flog,"Error while reading input file. Exiting.\n");
		}
		if (verbose_level>0) error_quit(" ** Error while reading input file. Exiting. **\n");
		else return 1;
	}
	if(procId == 0) {
		if (flog){
			fprintf(flog,"Fixed focal mechanism: [str, dip]=[%.3lf, %.3lf].\n", (*crst).str0[0],(*crst).dip0[0]);
			fflush(flog);
		}
	}

	//--------------read grid file:------------------------------//

	err = read_csep_template(fnametemplate, &no_magbins, &((*crst).nLat_out), &((*crst).nLon_out),
							 &((*crst).nD_out), &((*crst).N_allP), &((*crst).dlat_out), &((*crst).dlon_out),
							 &((*crst).ddepth_out), &((*crst).dmags), &olats, &olons, &odeps, 0,
							 &((*crst).latmin), &((*crst).latmax), &((*crst).lonmin), &((*crst).lonmax),
							 &((*crst).depmin), &((*crst).depmax), &mag1, &mag2, &((*crst).uniform));

	if(procId == 0) {
		if (flog && err) fprintf(flog, "Error while reading grid template (%s). Exiting.\n", fnametemplate);
	}
	if (verbose_level>0 && err!=0) error_quit(" ** Error while reading grid template. Exiting. **\n");

	lat0=(*crst).latmin;
	lat1=(*crst).latmax;
	lon0=(*crst).lonmin;
	lon1=(*crst).lonmax;
	d0=(*crst).depmin;
	d1=(*crst).depmax;

	(*crst).lat_out=olats;
	(*crst).lon_out=olons;
	(*crst).depth_out=odeps;
	(*crst).lat0=0.5*(lat1+lat0);
	(*crst).lon0=0.5*(lon1+lon0);

	if(procId == 0) {
		if (flog) {
			fprintf(flog, "Model domain: \n lat=[%.2lf, %.2lf], %d points; \n lon=[%.2lf, %.2lf], %d points; \n dep=[%.2lf, %.2lf], %d points; \n",
					lat0, lat1, (*crst).nLat_out, lon0, lon1, (*crst).nLon_out, d0, d1, (*crst).nD_out);
			fprintf(flog, " %s grid found.\n", ((*crst).uniform)? "Uniform" : "Non uniform");
		}
	}

	//--------------calculate magnitude bins:-------------//

	(*crst).nmags=no_magbins;
	(*crst).mags=dvector(1,no_magbins);

	for (int i=1; i<=no_magbins; i++) (*crst).mags[i]=mag1+(i-1)*(*crst).dmags;

	if(procId == 0) {
		if (flog) fprintf(flog, " mag=[%.2lf, %.2lf], %d bins.\n", (*crst).mags[1], (*crst).mags[(*crst).nmags], (*crst).nmags);
	}

	//--------------calculate refined geometry:-------------//

	dy=Re*DEG2RAD*(*crst).dlat_out;
	dx=Re*DEG2RAD*cos(DEG2RAD*(0.5*(lat1+lat0)))*(*crst).dlon_out;
	no_subpointsx= (int) (0.01+ceil(dx/resxy));
	no_subpointsy= (int) (0.01+ceil(dy/resxy));
	no_subpointsz= (int) (0.01+ceil((*crst).ddepth_out/resz));

	is_refined= (no_subpointsx*no_subpointsy*no_subpointsz>1) & (*crst).uniform;

	if ((*crst).uniform){
		NLat=(*crst).nLat_out;
		NLon=(*crst).nLon_out;
		Nd=(*crst).nD_out;

		NLat*=no_subpointsy;
		NLon*=no_subpointsx;
		Nd*=no_subpointsz;
		(*crst).nLat=NLat;
		(*crst).nLon=NLon;
		(*crst).nD=Nd;
		(*crst).N_allP=NG=(*crst).nLat*(*crst).nLon*(*crst).nD;

		//assume that lat0, lat1 are the boundaries of the domain (*not* the coordinates of the outermost cell centers).
		(*crst).dlat=(lat1-lat0)/(*crst).nLat;
		(*crst).dlon=(lon1-lon0)/(*crst).nLon;
		(*crst).ddepth=(d1-d0)/(*crst).nD;
		(*crst).lat=dvector(1,NG);
		(*crst).lon=dvector(1,NG);
		(*crst).depth=dvector(1,NG);
		(*crst).x=dvector(1,NG);
		(*crst).y=dvector(1,NG);
		(*crst).dAgrid=dvector(1,NG);
		(*crst).list_allP=ivector(1,NG);
		for (int i=1; i<=NG; i++) (*crst).list_allP[i]=i;

		for (int d=1; d<=Nd; d++){
			for (int lo=1; lo<=NLon; lo++){
				for (int la=1; la<=NLat; la++){
					ind=(d-1)*NLat*NLon+(lo-1)*NLat+la;
					(*crst).lat[ind]=lat0+(la-0.5)*(*crst).dlat;
					(*crst).lon[ind]=lon0+(lo-0.5)*(*crst).dlon;
					(*crst).depth[ind]=d0+(d-0.5)*(*crst).ddepth;
				}
			}
		}
	}
	else {
		if (no_subpointsx!=1 || no_subpointsy!=1 || no_subpointsz!=1 && verbose_level>0) {
			if(procId == 0) {
				if (verbose_level>1) printf("** Warning: non uniform grid in file %s, can not refine geometry (read_crust.c).**\n",fnametemplate);
				if (flog) fprintf(flog,"** Warning: non uniform grid in file %s, can not refine geometry (read_crust.c).**\n",fnametemplate);
			}
		}
		(*crst).nLat=0;
		(*crst).nLon=0;
		(*crst).nD=0;
		NG=(*crst).N_allP;

		//assume that lat0, lat1 are the boundaries of the domain (*not* the coordinates of the outermost cell centers).
		(*crst).dlat=(*crst).dlat_out;
		(*crst).dlon=(*crst).dlon_out;
		(*crst).ddepth=(*crst).ddepth_out;
		(*crst).lat=(*crst).lat_out;
		(*crst).lon=(*crst).lon_out;
		(*crst).depth=(*crst).depth_out;
		(*crst).x=dvector(1,NG);
		(*crst).y=dvector(1,NG);
		(*crst).dAgrid=dvector(1,NG);
		(*crst).list_allP=ivector(1,NG);
		for (int i=1; i<=NG; i++) (*crst).list_allP[i]=i;
	}

	if(procId == 0) {
		if (flog) {
			fprintf(flog, "Forecast resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km;\n", dy, dx, (*crst).ddepth_out);
			if (is_refined) {
				fprintf(flog, "Internal resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km -> %d x %d x %d = %d grid points.\n", resxy, resxy, resz, NLat, NLon, Nd, NG);
			}
			else {
				fprintf(flog, "Internal resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km-> %d x %d x %d = %d grid points.\n", dy, dx, (*crst).ddepth_out, (*crst).nLat_out, (*crst).nLon_out, (*crst).nD_out);
				fprintf(flog, "Real int.resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km.\n", dy, dx, (*crst).ddepth_out);
			}
		}
	}

	//--------------calculate area of each grid cell, and local coordinates:-------------//

	dAeq=pow(Re*PI/180,2)*(*crst).dlon*(*crst).dlat;
	for (int k=1; k<=NG;k++){
		(*crst).dAgrid[k]= dAeq*cos((*crst).lat[k]*PI/180);
		latlon2localcartesian((*crst).lat[k], (*crst).lon[k], (*crst).lat0, (*crst).lon0, (*crst).y+k, (*crst).x+k);
	}

	//--------------read value of focal mechanism grid from file:-------------//

	if (focmecgridfile && strcmp(focmecgridfile,"")!=0){

		err1 = read_focmecgridfile(focmecgridfile, crst);

		if (is_refined){
			err1+=convert_geometry((*crst),(*crst).str0, &strtmp, 0, 1);
			(*crst).str0=strtmp;
			err1+=convert_geometry((*crst),(*crst).dip0, &diptmp, 0, 1);
			(*crst).dip0=strtmp;
			err1+=convert_geometry((*crst),(*crst).rake0, &raketmp, 0, 1);
			(*crst).rake0=strtmp;
		}

		if(err1 && procId == 0) {
			printf("*Warning: errors occurred while reading focmecgridfile (%s)*\n", focmecgridfile);
			if (flog) {
				fprintf(flog, "*Warning: errors occurred while reading focmecgridfile*\n");
				fflush(flog);
			}
		}
	}

	if(procId == 0) {
		if (verbose_level>0)  printf("done\n");
	}

	return(err!=0);
}

int read_farfalle_crust(char * file, struct crust *crst){

	FILE *fin;
	int Nchar=200, err=0, junk;
	char line[Nchar];
	double s[3];	//regional stress field description;
	double st[3];	//regional stress field description;
	double di[3];	//regional stress field description;

	int procId = 0;
	int fileError = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	if(procId == 0) {
		fin=fopen(file,"r");
		if (!fin) {
			if (verbose_level) printf("**Error: can not find input file %s (read_farfalle_crust).**", file);
			if (flog) fprintf(flog, "**Error: can not find input file %s (read_farfalle_crust).**", file);

			fileError = 1;
		}

		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin); //useless field (controls Farfalle output).
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
		err+=ferror(fin);
		sscanf(line,"%lf %lf", &((*crst).lambda), &((*crst).mu));
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin); //useless field ('Adding regional stress field').
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
		err+=ferror(fin);
		sscanf(line,"%lf %lf %lf", s, s+1, s+2);
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
		err+=ferror(fin);
		sscanf(line,"%lf %lf", st, di);
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
		err+=ferror(fin);
		sscanf(line,"%lf %lf", st+1, di+1);
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
		err+=ferror(fin);
		sscanf(line,"%lf %lf", st+2, di+2);
		for (int i=1; i<=8; i++){
			line[0]='!';
			while (line[0]=='!') fgets(line,Nchar,fin);	//useless fields.
		}
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
		sscanf(line,"%lf %lf", &((*crst).fric), &((*crst).skepton));
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
	sscanf(line,"%d %lf", &junk, (*crst).str0);
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
	sscanf(line,"%d %lf", &junk, (*crst).dip0);
		line[0]='!';
		while (line[0]=='!') fgets(line,Nchar,fin);
	sscanf(line,"%d %lf", &junk, (*crst).rake0);

		fclose(fin);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&((*crst).lambda), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&((*crst).mu), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&((*crst).fric), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&((*crst).skepton), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&((*crst).str0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&((*crst).dip0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&((*crst).rake0), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(s, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(st, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(di, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	(*crst).S=prestress_eigen(s, st, di);

	return (err!=0);

}

int read_pscmp_crust(char *fname, struct crust *crst){
	FILE *fin;
	int dumerror=0;
	double junk;
	int nchar=200;
	char line[nchar];
	char comm[]="#";
	double s1, s2, s3;	//regional stress field description (see Wang input file).

	int procId = 0;
	int fileError = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	(*crst).lambda=31226, (*crst).mu=26624;//calculated for Vp=5.7,Vs=3.2, rho=2600 (from Wang psgrn input file for Parkfield). MPa.

	if(procId == 0) {
		if (verbose_level>1) printf("Loading model setup...");

		if (!(fin=fopen(fname,"r"))) {
			if (verbose_level) printf("Error: can not open file %s (read_pscmp_crust), Exiting.\n", fname);
			if (flog) fprintf(flog, "Error: can not open file %s (read_pscmp_crust), Exiting.\n", fname);

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
		line[0]=comm[0];
		while(line[0]==comm[0]) fgets(line,nchar,fin);
		while(line[0]!=comm[0]) fgets(line,nchar,fin);
		while(line[0]==comm[0]) fgets(line,nchar,fin);
		fgets(line,nchar,fin);
		dumerror = sscanf(line, " %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf",
						  &junk, &((*crst).fric), &((*crst).skepton),
						  (*crst).str0, (*crst).dip0, (*crst).rake0,
						  &s1, &s2, &s3);
		fclose(fin);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&dumerror, 1, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Bcast(&((*crst).fric), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&((*crst).skepton), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*crst).str0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*crst).dip0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast((*crst).rake0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&s1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&s2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&s3, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	prestress(s1, s2, s3, (*crst).str0[0], (*crst).dip0[0], (*crst).rake0[0], 0.0,(*crst).fric, &((*crst).S));

	return (dumerror!=9);
}

int read_focmecgridfile(char *fname, struct crust *crst) {
	// Variables used for MPI
	int procId = 0;
	int fileError = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double strtmp, diptmp, raketmp;
	double **data;
	int err=0, NL, NP=(*crst).nLat_out*(*crst).nLon_out*(*crst).nD_out;

	data=dmatrix(1,2,1,(*crst).N_allP);

	if(procId == 0) {
		err = read_matrix(fname, 2, 0, data, &NL);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NL, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&err, 1, MPI_LONG, 0, MPI_COMM_WORLD);

		long nrl=1, nrh=2, ncl=1, nch=(*crst).N_allP;
		long nrow=nrh-nrl+1, ncol=nch-ncl+1;
		MPI_Bcast(data[nrl], nrow*ncol+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	if (err || NL!=NP){
		if (NL!=NP) {
			if(procId == 0) {
				if (verbose_level) printf("Error: wrong number of lines in file %s (%d lines found; %d expected). (read_focmecgridfile).\n", fname, NL, NP);
				if (flog) fprintf(flog, "Error: wrong number of lines in file %s (%d lines found; %d expected). (read_focmecgridfile).\n", fname, NL, NP);
			}
		}
		else {
			if(procId == 0) {
				if (verbose_level) printf("Error: can not open file %s (read_focmecgridfile), exiting.\n", fname);
				if (flog) fprintf(flog, "Error: can not open file %s (read_focmecgridfile), exiting.\n", fname);
			}
		}

		free_dmatrix(data, 1,2,1,(*crst).N_allP);

		return 1;
	}

	strtmp=(*crst).str0[0];
	diptmp=(*crst).dip0[0];
	//raketmp=(*crst).rake0[0];
	(*crst).variable_fixmec=1;
	(*crst).str0=dvector(0,(*crst).N_allP);
	(*crst).dip0=dvector(0,(*crst).N_allP);
	//(*crst)->rake0=dvector(0,crst.N_allP);

	(*crst).str0[0]=strtmp;
	(*crst).dip0[0]=diptmp;
	for (int i=1; i<=crst->N_allP; i++){
		(*crst).str0[i]=data[1][i];
		(*crst).dip0[i]=data[2][i];
	}

	free_dmatrix(data, 1,2,1,(*crst).N_allP);

	return 0;
}

