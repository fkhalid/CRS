#include "print_output.h"

//void printCSEPforecast(char *filename, double *lats, double *lons, double *deps, double *mags, double *rates, double *mag_fact, int NG, int Nmag){
//
//	double Lat1,Lat2,Lon1,Lon2, D1, D2;
//	double dlat, dlon, ddep;
//	FILE *fout;
//	fout=fopen(filename,"w");
//	counter=1;
//	while (lons[counter]==lons[1]) counter++;
//	dlon= lons[counter]-lons[1];
//	counter=1;
//	while (lats[counter]==lats[1]) counter++;
//	dlat= lats[counter]-lats[1];
//	counter=1;
//	while (deps[counter]==deps[1]) counter++;
//	ddep= deps[counter]-deps[1];
//
//	for (int m=1; m<=Nmag; m++){
//		for (int k=1; k<=NG; k++){
//			Lon1=lons[k]-dlon/2.0;
//			Lon2=lons[k]+dlon/2.0;
//			Lat1=lats[k]-dlat/2.0;
//			Lat2=lats[k]+dlat/2.0;
//			D1=depths[k]-ddep/2.0;
//			D2=depths[k]+ddep/2.0;
//			M1=mags[m]-dmag/2.0;
//			M2=mags[m]+dmag/2.0;
//			fprintf(foutLoop2,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, M1, M2, rates[k]*mag_fact[m]);
//		}
//	}
//	fclose(fout);
//
//	return;
//}

int sum_DCFS(struct pscmp *DCFS, double **cmb, int N, int Ntot){
//adds up all DCFS[0...N-1], and returns vector cmb containing cumulative field.
//Ntot: total no. of grid points.
//if *cmb==NULL memory is allocated, otherwise not (make sure vector has correct no. of elements).

	int i;

	if (*cmb==NULL) *cmb=dvector(1,Ntot);
	for (int k=1; k<=Ntot; k++) (*cmb)[k]=0.0;


	for (int n=0; n<N; n++){
		i=1;
		for (int k=1; k<=Ntot; k++){		//todo check if this works with DCFS.nsel<Ntot.
			if (i<=DCFS[n].nsel && DCFS[n].which_pts[i]==k) {
				(*cmb)[k]+=DCFS[n].cmb[i];
				i++;
			}
		}
	}

	return 0;

}

void printmatrix(double **S){
	//printout 3x3 matrix

	for (int ii=1; ii<=3; ii++) {
		for (int jj=1; jj<3; jj++) printf("%.3lf\t",S[ii][jj]);
		if (ii<3) printf("%.3lf\n",S[ii][3]);
		else  printf("%.3lf]\n",S[ii][3]);
	}

	return;

}

int print_cat(char *fname, struct catalog cat){
/*uses following format:
 * 1. Longitude [deg]
   2. Latitude [deg]
   3. Time (days)
   4. Magnitude
   5. Depth [km]
 *
 */

	FILE *fout;
	int Z=cat.Z;

	fout=fopen(fname,"w");
	for (int i=1; i<Z; i++) fprintf(fout, "%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf%.3lf\t%.3lf\t\n", cat.lon0[i], cat.lat0[i], cat.t[i], cat.mag[i], cat.depths0[i], cat.err[i], cat.verr[i]);
	fprintf(fout, "%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf%.3lf\t%.3lf\t", cat.lon0[Z], cat.lat0[Z], cat.t[Z], cat.mag[Z], cat.depths0[Z],cat.err[Z], cat.verr[Z]);
	fclose(fout);

	return 0;
}

int print_cat_long(char *fname, struct catalog cat){
/*uses following format (2 lines per event):
 * No_gridpoints herr verr point1 point2 point3 ...
 * No_gridpoints herr verr weight1 weight2 weight3 ...
 */

	FILE *fout;
	int Z=cat.Z;

	fout=fopen(fname,"w");
	for (int i=1; i<Z; i++) {

		fprintf(fout, "%d\t%.3e\t%.3e\t", cat.ngrid[i], cat.err[i], cat.verr[i]);
		for (int j=1; j<=cat.ngrid[i]; j++) fprintf(fout, "%d\t", cat.ngridpoints[i][j]);

		fprintf(fout, "\n%d\t%.3e\t%.3e\t", cat.ngrid[i], cat.err[i], cat.verr[i]);
		for (int j=1; j<=cat.ngrid[i]; j++) fprintf(fout, "%.3e\t", cat.weights[i][j]);
		fprintf(fout, "\n");
	}

	fclose(fout);

	return 0;
}

int print_rate(char *fname, struct crust crst, double Mmin, double *rate){
/* line print_grid, but doesn't use struct pscmp.
 * if rate==NULL, will printout (crst.r0)*(crst.rate0)/Ntot.
 */

		double Lon1, Lon2, Lat1, Lat2, D1, D2;
		FILE *fout;
		int i, Ntot= crst.N_allP;
		double *r, r0;

		if (rate) {
			r=rate;
			r0=1.0*Ntot;	//since will printout rate/Ntot (due to definition of crst.rate0).
		}

		else {
			r=crst.rate0;
			r0=crst.r0;
		}


		fout=fopen(fname,"w");
		if (fout==NULL){
			if (verbose_level>0) printf("Error: file %s could not be opened (print_cmb). \n",fname);
			return(1);
		}
		i=1;
		for (int k=1; k<=Ntot; k++){		//todo check if this works with DCFS.nsel<Ntot.
			Lon1=crst.lon[k]-crst.dlon/2.0;
			Lon2=crst.lon[k]+crst.dlon/2.0;
			Lat1=crst.lat[k]-crst.dlat/2.0;
			Lat2=crst.lat[k]+crst.dlat/2.0;
			D1=crst.depth[k]-crst.ddepth/2.0;
			D2=crst.depth[k]+crst.ddepth/2.0;
			fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, Mmin, 8.0, r0*r[k]/(1.0*Ntot));
		}
		fclose(fout);

		return(0);
}


int print_grid(char *fname, struct pscmp DCFS, struct crust crst, double *rate){
/* if rate is a null pointer, prints out the coulomb stress field (DCFS.cmb).
 * uses refined grid geometry.
 */

	double *r;
	double Lon1, Lon2, Lat1, Lat2, D1, D2;
	FILE *fout;
	int i, Ntot=crst.uniform? (crst.nLat*crst.nLon*crst.nD) : crst.N_allP;

	r=(rate==NULL) ? DCFS.cmb : rate;

	fout=fopen(fname,"w");
	if (fout==NULL){
		if (verbose_level>0) printf("Error: file %s could not be opened (print_cmb). \n",fname);
		return(1);
	}
	i=1;
	for (int k=1; k<=Ntot; k++){		//todo check if this works with DCFS.nsel<Ntot.
		Lon1=crst.lon[k]-crst.dlon/2.0;
		Lon2=crst.lon[k]+crst.dlon/2.0;
		Lat1=crst.lat[k]-crst.dlat/2.0;
		Lat2=crst.lat[k]+crst.dlat/2.0;
		D1=crst.depth[k]-crst.ddepth/2.0;
		D2=crst.depth[k]+crst.ddepth/2.0;
		if (i<=DCFS.nsel && DCFS.which_pts[i]==k) {
			//rate would contain all elements, but DCFS.cmb only those selected.
			if (rate) fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, 3.0, 8.0, r[k]);
			else fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, 3.0, 8.0, r[i]);
			i++;
		}
		else fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, 3.0, 8.0, 0.0);
	}
	fclose(fout);

	return(0);
}

int print_slipmodel(char* filename, struct eqkfm *eqfm1, int NF){
// filename with extension; NF=no. of faults (elements of eqkfm).

	FILE *fout;
	int err=0;

	fout=fopen(filename,"w");
	if (fout==NULL) {
		if (verbose_level>1) printf("Warning: file %s could not be written (by function: print_slipmode).\n",filename);
		return 1;
	}
	else{
		for (int f=0; f<NF; f++){
			switch (eqfm1[f].whichfm){
				case 1:
					fprintf(fout, "%d   %.4lf   %.4lf   %.3lf   %.2lf   %.2lf   %.3lf   %.3lf   %d   %d   %.5lf\n", f+1, eqfm1[f].lat, eqfm1[f].lon, eqfm1[f].depth, eqfm1[f].L, eqfm1[f].W, eqfm1[f].str1, eqfm1[f].dip1, eqfm1[f].np_st, eqfm1[f].np_di, eqfm1[f].t);
					break;
				case 2:
					fprintf(fout, "%d   %.4lf   %.4lf   %.3lf   %.2lf   %.2lf   %.3lf   %.3lf   %d   %d   %.5lf\n", f+1, eqfm1[f].lat, eqfm1[f].lon, eqfm1[f].depth, eqfm1[f].L, eqfm1[f].W, eqfm1[f].str2, eqfm1[f].dip2, eqfm1[f].np_st, eqfm1[f].np_di, eqfm1[f].t);
					break;
				case 0:
					if (verbose_level>1) printf("Warning: whichfm=0, using first plane (print_slipmodel).\n");
					fprintf(fout, "%d   %.4lf   %.4lf   %.3lf   %.2lf   %.2lf   %.3lf   %.3lf   %d   %d   %.5lf\n", f+1, eqfm1[f].lat, eqfm1[f].lon, eqfm1[f].depth, eqfm1[f].L, eqfm1[f].W, eqfm1[f].str1, eqfm1[f].dip1, eqfm1[f].np_st, eqfm1[f].np_di, eqfm1[f].t);
					break;
				default:
					if (verbose_level>1) printf("Warning: ambiguos focal plane  -> output not written (print_slipmodel).\n");
					err+=1;
					continue;		//associated with for loop (not switch).
			}
			for (int p=1; p<=eqfm1[f].np_st*eqfm1[f].np_di; p++) {
				fprintf(fout, "%12.5lf\t%12.5lf\t%12.5lf\t%12.5lf\t%12.5lf\n", eqfm1[f].pos_s[p], eqfm1[f].pos_d[p], eqfm1[f].slip_str[p], eqfm1[f].slip_dip[p], 0.0);
			}
		}
		fclose(fout);
	}
	return err;
}
