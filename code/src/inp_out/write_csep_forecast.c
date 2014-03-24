#include "write_csep_forecast.h"

void csep_forecast(char *filename, struct crust crst, double *rates, int calculation_grid){

	double *lats, *lons, *deps, *mags, *magGR;
	double dlat, dlon, ddep, dmag;
	int NG, Nmag;

	if (calculation_grid){
		lats=crst.lat;
		lons=crst.lon;
		deps=crst.depth;
		dlat=crst.dlat;
		dlon=crst.dlon;
		ddep=crst.ddepth;
		NG=crst.uniform? (crst.nD*crst.nLat*crst.nLon) : crst.N_allP;
	}

	else{
	lats=crst.lat_out;
	lons=crst.lon_out;
	deps=crst.depth_out;
	dlat=crst.dlat_out;
	dlon=crst.dlon_out;
	ddep=crst.ddepth_out;
	NG=crst.uniform? (crst.nD_out*crst.nLat_out*crst.nLon_out) : crst.N_allP;
	}


	mags=crst.mags;
	magGR=crst.GRmags;
	dmag=crst.dmags;
	Nmag=crst.nmags;

	write_csep_forecast(filename, lats, lons, deps, dlat, dlon, ddep, mags, dmag, rates, magGR, NG, Nmag);
	return;
}

void write_csep_forecast(char *filename, double *lats, double *lons, double *deps, double dlat, double dlon, double ddep,
		double *mags, double dmag, double *rates, double *mag_fact, int NG, int Nmag){

	//input:
	//lats, lons, deps [1...NG]: central point of cells.
	//dlat, dlon, ddep, dmag: cell sizes.
	//rates [1...NG]: rate at each cell (sum over all magnitude bins).
	//mags [1...Nmag]: center of magnitude bins.
	//mag_fact [1...Nmag]: factor by which rate should be scaled in each magnitude bin (they should add up to 1).


	double Lat1,Lat2,Lon1,Lon2, D1, D2, M1, M2;

	FILE *fout;
	fout=fopen(filename,"w");

	for (int k=1; k<=NG; k++){
		for (int m=1; m<=Nmag; m++){
			Lon1=lons[k]-dlon/2.0;
			Lon2=lons[k]+dlon/2.0;
			Lat1=lats[k]-dlat/2.0;
			Lat2=lats[k]+dlat/2.0;
			D1=deps[k]-ddep/2.0;
			D2=deps[k]+ddep/2.0;
			M1=mags[m]-dmag/2.0;
			M2=mags[m]+dmag/2.0;
			fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, M1, M2, rates[k]*mag_fact[m]);
		}
	}
	fclose(fout);

	return;
}
