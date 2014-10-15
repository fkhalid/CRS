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

int read_crust(char *fnametemplate, char *focmecgridfile, struct crust *crst, double resxy, double resz, int multiple_focmecfiles){
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
	double *dumrate=NULL, *bg_rate=NULL;

	print_screen("Loading model setup...");
	print_logfile("\nEntering read_crust...\n");

	double *olats, *olons, *odeps;

	//--------------read grid file:------------------------------//

	err = read_csep_template(fnametemplate, &no_magbins, &((*crst).nLat_out), &((*crst).nLon_out),
							 &((*crst).nD_out), &((*crst).N_allP), &((*crst).dlat_out), &((*crst).dlon_out),
							 &((*crst).ddepth_out), &((*crst).dmags), &olats, &olons, &odeps, (multiple_focmecfiles) ? &dumrate : NULL,
							 &((*crst).latmin), &((*crst).latmax), &((*crst).lonmin), &((*crst).lonmax),
							 &((*crst).depmin), &((*crst).depmax), &mag1, &mag2, &((*crst).uniform));

	if(err) {
		print_logfile("Error while reading grid template (%s). Exiting.\n", fnametemplate);
		error_quit(" ** Error while reading grid template. Exiting. **\n");
	}

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

	print_logfile( "Model domain: \n lat=[%.2lf, %.2lf], %d points; \n lon=[%.2lf, %.2lf], %d points; \n dep=[%.2lf, %.2lf], %d points; \n",
					lat0, lat1, (*crst).nLat_out, lon0, lon1, (*crst).nLon_out, d0, d1, (*crst).nD_out);
	print_logfile(" %s grid found.\n", ((*crst).uniform)? "Uniform" : "Non uniform");

	//--------------calculate magnitude bins:-------------//

	(*crst).nmags=no_magbins;
	(*crst).mags=dvector(1,no_magbins);

	for (int i=1; i<=no_magbins; i++) (*crst).mags[i]=mag1+(i-1)*(*crst).dmags;

	print_logfile(" mag=[%.2lf, %.2lf], %d bins.\n", (*crst).mags[1], (*crst).mags[(*crst).nmags], (*crst).nmags);

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
		if (no_subpointsx!=1 || no_subpointsy!=1 || no_subpointsz!=1) {
			print_screen("** Warning: non uniform grid in file %s, can not refine geometry (read_crust.c).**\n",fnametemplate);
			print_logfile("** Warning: non uniform grid in file %s, can not refine geometry (read_crust.c).**\n",fnametemplate);
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

	print_logfile("Forecast resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km;\n", dy, dx, (*crst).ddepth_out);
	if (is_refined) {
		print_logfile( "Internal resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km -> %d x %d x %d = %d grid points.\n", resxy, resxy, resz, NLat, NLon, Nd, NG);
	}
	else {
		print_logfile("Internal resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km-> %d x %d x %d = %d grid points.\n", dy, dx, (*crst).ddepth_out, (*crst).nLat_out, (*crst).nLon_out, (*crst).nD_out, (*crst).N_allP);
		print_logfile("Real int.resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km.\n", dy, dx, (*crst).ddepth_out);
	}

	//--------------calculate area of each grid cell, and local coordinates:-------------//

	dAeq=pow(Re*PI/180,2)*(*crst).dlon*(*crst).dlat;
	for (int k=1; k<=NG;k++){
		(*crst).dAgrid[k]= dAeq*cos((*crst).lat[k]*PI/180);
		latlon2localcartesian((*crst).lat[k], (*crst).lon[k], (*crst).lat0, (*crst).lon0, (*crst).y+k, (*crst).x+k);
	}

	//--------------read value of focal mechanism grid from file:-------------//
	//----------(this applies is a single foc mec is given per grid point)----//

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

		if(err1){
			print_screen("*Warning: errors occurred while reading focmecgridfile (%s)*\n", focmecgridfile);
			print_logfile("*Warning: errors occurred while reading focmecgridfile*\n");
		}
	}

	//---------------read indices and no. of foc. mec zones into crst-----------//
	//---(this applies is a set of foc mec is associated to each grid point)----//

	else{
		if(multiple_focmecfiles) {
			err+=convert_geometry(*crst, dumrate, &bg_rate, 1, 1);
			for (int n=1; n<=(*crst).N_allP; n++) (*crst).fmzone[n]= (int)bg_rate[n] -1;
			(*crst).nofmzones=(int) max_v(bg_rate+1,(*crst).N_allP);
			free_dvector(dumrate, 1, 1);
			free_dvector(bg_rate, 1, 1);
		}
	}


	print_screen("done\n");

	return(err!=0);
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
			print_screen("Error: wrong number of lines in file %s (%d lines found; %d expected). (read_focmecgridfile).\n", fname, NL, NP);
			print_logfile("Error: wrong number of lines in file %s (%d lines found; %d expected). (read_focmecgridfile).\n", fname, NL, NP);
		}
		else {
			print_screen("Error: can not open file %s (read_focmecgridfile), exiting.\n", fname);
			print_logfile("Error: can not open file %s (read_focmecgridfile), exiting.\n", fname);
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

