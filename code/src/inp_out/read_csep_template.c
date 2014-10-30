/*
 * read_txttemplate.c
 *
 *  Created on: Oct 4, 2013
 *      Author: camcat
 */

#include "read_csep_template.h"

#include <math.h>
#include <stddef.h>

#include "../defines.h"
#include "../geom/convert_geometry.h"
#include "../util/nrutil.h"
#include "read_matrix.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int read_rate(struct crust crst, char *fname, double **bg_rate, double *minmag){
	/* Read 9th column of a forecast-like file (same structure as output grid file), and converts it to internal geometry.
	 *
	 * Input:
	 * 	fname: file name with structure of output file (or forecast template)
	 * 	crst: crst structure used to convert geometry
	 *
	 * Ouput:
	 * 	bg_rate: rate (col. 9 of file), reshaped with internal geometry.
	 * 	minmag: min. mag bin found in file.
	 */

	double *dum_rate, dmag;
	int err;

	err=read_csep_template(fname, 0,0,0,0,0,0,0,0,&dmag,0,0,0,&dum_rate,0,0,0,0,0,0,minmag,0, NULL);
	err+=convert_geometry(crst, dum_rate, bg_rate, 1, 1);
	if (minmag) *minmag-=0.5*dmag;
	return err;
}

int read_csep_template(char *fname, int *no_magbins, int *nlat, int *nlon,
					   int *ndep, int *ng, double *dlat, double *dlon,
					   double *ddep, double *dmag, double **lats, double **lons,
					   double **deps, double **rate, double *minlat, double *maxlat,
					   double *minlon, double *maxlon, double *mindep, double *maxdep,
					   double *minmag, double *maxmag, int *uni) {

/* Reads a template file in csep format.
 *
 * Input:
 * 	fname, name of txt file.
 *
 * Output:
 * 	no_magbins: no. of magnitude bins;
 * 	nlat, nlon ndep: no of gridpoints with different lat, lon, dep (may not work for non homogeneous grid).
 * 	ng: tot no of grid points.
 * 	dlat, dlon, ddep, dmag: shortst distance between grid points (for lat, lon, this is the "defaultCellDimension" value; for depth, it is calculated).
 * 	lats, lons, deps: 1D arrays containing grid points coordinates: [1...ng].
 * 	minX, maxX: edges of domain (lat, lon, depth);
 * 	minmag, maxmag: smallest, largest *centers* of magnitude bin.
 * 	rate: rate in each cell (summed over magnitude bins).
 *
 * All output variables can be passed as null pointers, and will be ignored.
 *
 * */

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int NC, NL, NH=0, NP;
	int n1, nmag, err=0;
	long NR;
	double **data;
	double lat0=1e30, lat1=-1e30, lon0=1e30, lon1=-1e30, dep0=1e30, dep1=-1e30, mag0, mag1;
	double dlati, dloni, ddepi, dmagi;
	double toll=1e-6;
	double closest_lat, closest_lon, closest_dep;

	if(procId == 0) {
		NC = countcol(fname);
		NL = countline(fname);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NC, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NL, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (NL>0 && NC>0) {
		data = dmatrix(1,NC, 1, NL+1);

		if(procId == 0) {
			print_screen("Reading template file %s\n", fname);
			print_logfile("Reading template file %s\n", fname);
			err = read_matrix(fname, NC, NH, data, &NR);
		}

		#ifdef _CRS_MPI
			MPI_Bcast(&err, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		#endif

		if(err) {
			return 1;
		}

		#ifdef _CRS_MPI
			MPI_Bcast(&NR, 1, MPI_LONG, 0, MPI_COMM_WORLD);

			long nrl=1, nrh=NC, ncl=1, nch=NL+1;
			long nrow=nrh-nrl+1, ncol=nch-ncl+1;
			MPI_Bcast(data[nrl], nrow*ncol+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif
	}
	else {
		return 1;
	}

	nmag=2;
	mag0=data[7][1];
	mag1=data[8][1];
	// todo [coverage] this block is never tested
	while (fabs(data[8][1]-data[8][nmag])>toll) {
		mag0=fmin(mag0, data[7][nmag]);
		mag1=fmax(mag0, data[8][nmag]);
		nmag++;
	}
	nmag-=1;

	NP=NR/nmag;
	//assume cells have all the same size.
	dlati=data[4][1]-data[3][1];
	dloni=data[2][1]-data[1][1];
	ddepi=data[6][1]-data[5][1];
	dmagi=data[8][1]-data[7][1];

	if (lats) *lats=dvector(1,NP);
	if (lons) *lons=dvector(1,NP);
	if (deps) *deps=dvector(1,NP);
	if (rate) *rate=dvector(1,NP);

	closest_lat=1e30;
	closest_lon=1e30;
	closest_dep=1e30;

	for (int n=1; n<=NP; n++){
		n1=nmag*n;
		if (lats) (*lats)[n]= 0.5*(data[3][n1]+data[4][n1]);
		if (lons) (*lons)[n]= 0.5*(data[1][n1]+data[2][n1]);
		if (deps) (*deps)[n]= 0.5*(data[5][n1]+data[6][n1]);
		if (rate) {
			(*rate)[n]=0.0;
			for (int i=n1-nmag+1; i<=n1; i++)(*rate)[n]+= data[9][i];
		}

		lat0=fmin(lat0, data[3][n1]);
		lat1=fmax(lat1, data[4][n1]);
		lon0=fmin(lon0, data[1][n1]);
		lon1=fmax(lon1, data[2][n1]);
		dep0=fmin(dep0, data[5][n1]);
		dep1=fmax(dep1, data[6][n1]);

		if (fabs(data[3][n1]-data[3][1])>toll) closest_lat=fmin(closest_lat,fabs(data[3][n1]-data[3][1]));
		if (fabs(data[1][n1]-data[1][1])>toll) closest_lon=fmin(closest_lon,fabs(data[1][n1]-data[1][1]));
		if (fabs(data[5][n1]-data[5][1])>toll) closest_dep=fmin(closest_dep,fabs(data[5][n1]-data[5][1]));

	}

	//calculated for uniform grid: (toll since casting is same as floor).
	//fixme testing for uniform grid gives false positives!
	if (fabs(closest_lat-dlati)<toll && fabs(closest_lon-dloni)<toll && fabs(closest_dep-ddepi)<toll) {

		if (uni)  *uni=1;
		if (nlat) *nlat= (int) (toll+(lat1-lat0)/dlati);
		if (nlon) *nlon= (int) (toll+(lon1-lon0)/dloni);
		if (ndep) *ndep= (int) (toll+(dep1-dep0)/ddepi);
	}

	// todo [coverage] this block is never tested
	else {
		if (uni)  *uni=0;
		if (nlat) *nlat= 0;
		if (nlon) *nlon= 0;
		if (ndep) *ndep= 0;
	}

	if (no_magbins) *no_magbins = nmag;
	if (ng) *ng= NP;

	if (dlat) *dlat=dlati;
	if (dlon) *dlon=dloni;
	if (ddep) *ddep=ddepi;
	if (dmag) *dmag=dmagi;

	if (minlat) *minlat=lat0;
	if (maxlat) *maxlat=lat1;
	if (minlon) *minlon=lon0;
	if (maxlon) *maxlon=lon1;
	if (mindep) *mindep=dep0;
	if (maxdep) *maxdep=dep1;
	if (minmag) *minmag=mag0+0.5*dmagi;
	if (maxmag) *maxmag=mag1-0.5*dmagi;

	if(procId == 0) {
		NC = countcol(fname);
		NL = countline(fname);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NC, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NL, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	free_dmatrix(data,1,NC, 1, NL+1);

	return (0);
}


//int read_csep_template_uni(char *fname, int *no_magbins, int *nlat, int *nlon, int *ndep, int *ng, double *dlat, double *dlon, double *ddep, double *dmag,
//		double **lats, double **lons, double **deps, double **rate, double *minlat, double *maxlat, double *minlon, double *maxlon, double *mindep, double *maxdep,
//		double *minmag, double *maxmag){
//
///* input:
// * fname, name of txt file.
// *
// * output:
// * no_magbins: no. of magnitude bins;
// * nlat, nlon ndep: no of gridpoints with different lat, lon, dep (may not work for non homogeneous grid).
// * ng: tot no of grid points.
// * dlat, dlon, ddep, dmag: shortst distance between grid points (for lat, lon, this is the "defaultCellDimension" value; for depth, it is calculated).
// * lats, lons, deps: 1D arrays containing grid points coordinates: [1...ng].
// * minX, maxX: edges of domain (lat, lon, depth);
// * minmag, maxmag: smallest, largest *centers* of magnitude bin.
// * rate: rate in each cell (summed over magnitude bins).
// *
// * All output variables can be passed as null pointers, and will be ignored.
// *
// * */
//
//	int NC=countcol(fname), NL=countline(fname), NH=0, NP;
//	int n1, nmag;
//	long NR;
//	double **data;
//	double lat0=1e30, lat1=-1e30, lon0=1e30, lon1=-1e30, dep0=1e30, dep1=-1e30, mag0, mag1;
//	double dlati, dloni, ddepi, dmagi;
//	double toll=1e-6;
//
//	if (NL>0 && NC>0) data=dmatrix(1,NC, 1, NL+1);
//	else return 1;
//	read_matrix(fname, NC, NH, data, &NR);
//
//	nmag=2;
//	mag0=data[7][1];
//	mag1=data[8][1];
//	while (fabs(data[8][1]-data[8][nmag])>toll) {
//		mag0=fmin(mag0, data[7][nmag]);
//		mag1=fmax(mag0, data[8][nmag]);
//		nmag++;
//	}
//	nmag-=1;
//
//	NP=NR/nmag;
//	//assume cells have all the same size.
//	dlati=data[4][1]-data[3][1];
//	dloni=data[2][1]-data[1][1];
//	ddepi=data[6][1]-data[5][1];
//	dmagi=data[8][1]-data[7][1];
//
//	if (lats) *lats=dvector(1,NP);
//	if (lons) *lons=dvector(1,NP);
//	if (deps) *deps=dvector(1,NP);
//	if (rate) *rate=dvector(1,NP);
//
//	for (int n=1; n<=NP; n++){
//		n1=nmag*n;
//		if (lats) (*lats)[n]= 0.5*(data[3][n1]+data[4][n1]);
//		if (lons) (*lons)[n]= 0.5*(data[1][n1]+data[2][n1]);
//		if (deps) (*deps)[n]= 0.5*(data[5][n1]+data[6][n1]);
//		if (rate) {
//			(*rate)[n]=0.0;
//			for (int i=n1-nmag+1; i<=n1; i++)(*rate)[n]+= data[9][i];
//		}
//
//		lat0=fmin(lat0, data[3][n1]);
//		lat1=fmax(lat1, data[4][n1]);
//		lon0=fmin(lon0, data[1][n1]);
//		lon1=fmax(lon1, data[2][n1]);
//		dep0=fmin(dep0, data[5][n1]);
//		dep1=fmax(dep1, data[6][n1]);
//	}
//
//	//calculated for uniform grid: (toll since casting is same as floor).
//	if (nlat) *nlat= (int) (toll+(lat1-lat0)/dlati);
//	if (nlon) *nlon= (int) (toll+(lon1-lon0)/dloni);
//	if (ndep) *ndep= (int) (toll+(dep1-dep0)/ddepi);
//	if (no_magbins) *no_magbins = nmag;
//	if (ng) *ng= NP;
//
//	if (dlat) *dlat=dlati;
//	if (dlon) *dlon=dloni;
//	if (ddep) *ddep=ddepi;
//	if (dmag) *dmag=dmagi;
//
//	if (minlat) *minlat=lat0;
//	if (maxlat) *maxlat=lat1;
//	if (minlon) *minlon=lon0;
//	if (maxlon) *maxlon=lon1;
//	if (mindep) *mindep=dep0;
//	if (maxdep) *maxdep=dep1;
//	if (minmag) *minmag=mag0+0.5*dmagi;
//	if (maxmag) *maxmag=mag1-0.5*dmagi;
//
//	free_dmatrix(data,1,countcol(fname), 1, countline(fname)+1);
//
//	return (0);
//}
