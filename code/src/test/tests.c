/*
 * tests.c
 *
 *  Created on: Aug 23, 2013
 *      Author: camcat
 */


#include "tests.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//#include <unistd.h>

#include "../defines.h"
//#include "../general/forecast_stepG.h"
#include "../general/mem_mgmt.h"
#include "../geom/convert_geometry.h"
#include "../geom/coord_trafos.h"
#include "../geom/dist2fault.h"
#include "../geom/find_gridpoints.h"
#include "../inp_out/print_output.h"
#include "../inp_out/read_crust.h"
#include "../inp_out/read_csep_template.h"
#include "../inp_out/read_eqkfm.h"
//#include "../inp_out/read_eqkfm_fsp.h"
#include "../inp_out/read_focmec.h"
#include "../inp_out/read_inputfile.h"
#include "../inp_out/read_matrix.h"
#include "../inp_out/read_zmap.h"
#include "../inp_out/write_csep_forecast.h"
#include "../okada/okadaDCFS.h"
#include "../okada/prestress.h"
#include "../seis/background_rate.h"
#include "../seis/cmbopt.h"
#include "../seis/decluster.h"
#include "../seis/GR.h"
#include "../seis/Helmstetter.h"
#include "../seis/soumod1.h"
#include "../util/hash.h"
#include "../util/moreutil.h"
#include "../util/nr.h"
#include "../util/nrutil.h"

char *testfolder="test";

/* Tests to do:
 *
 * * choose which_taper in multiple fault systems with vertical alignment.
 * * write sample slip model files and run tests_eqkfm_addslipmodels().
 * * should change taper condition so that taper is not applied if dips are very different!
 * * all variants of okadaDCFS should give same output (when used in right order) --> check.
 *
 */

//int setup_catalogetc(char *catname, char *focmeccat, int NgridT, double t0, double tw1, double tw2, double t1, double dDCFS, double Mag, double Mag_main, struct pscmp DCFS0,
//		struct crust crst, struct catalog *cat, struct eqkfm **eqkfm1, int aftershocks, int *Nfm, int *Ntot, int *Nmain){}

//extern int xmlLoadExtDtdDefaultValue;
extern int gridPMax;
double DCFS_cap=1e7;

int test_readmultiplefocmec(){

	int nf=3;
	char **files;
	int NFM, NFM2;
	double **focmec;
	struct eqkfm *eqkfm;
	int NFM0, NFM20;
	double **focmec0;
	struct eqkfm *eqkfm0;
	struct crust crst;
	struct tm reftime;
	char *crust_file="/home/des/camcat/Code/CRS_2.01/input/Tohoku_new.inp";
	char *temp_file="/home/des/camcat/Code/CRS_2.01/input/tohoku_template_nonuniform2.dat";
	double res=1000.0;
	int *firstel;

	files=malloc(nf*sizeof(char *));
	for (int i=1; i<=nf; i++){
		files[i-1]=malloc(120*sizeof(char));
		sprintf(files[i-1],"input/other/NIED_GMTformat%d.dat",i);
	}

	sscanf("2011-03-11T14:46:18Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	reftime.tm_year-=1900;
	reftime.tm_mon-=1;
	reftime.tm_isdst=0;

	//broken (no crust_file anymore)
//	read_crust(crust_file, temp_file, &crst, NULL, res, res);
	readmultiplefocmec(files, nf, crst, 0.0, 0.0, 0.0, reftime, 0.0, 100, 100, 2.0, &focmec, &firstel, &NFM, &NFM2, &eqkfm, 1, 0);
	readfocmec("input/other/NIED_GMTformat.dat", crst, 0.0, 0.0, 0.0, reftime, 0.0, 100, 100, 2.0, &focmec0, &NFM0, &NFM20, &eqkfm0, 1, 0);

	for (int i=0; i<=nf; i++) printf("%d\n", firstel[i]);

	return 0;
}

int test_next_separator(){

	char fname[]="/home/des/camcat/Data/slipmodels/Tohoku/srcmod/s2011TOHOKU01IDEx.fsp";
	int nf;
	struct eqkfm *eqfm;

	read_fsp_eqkfm(fname, &eqfm, &nf);

	print_slipmodel("/home/des/camcat/Data/slipmodels/Tohoku/srcmod/s2011TOHOKU01IDEx", eqfm, nf);

//	adjust_faults(eqfm, nf, 0);
//
//	print_slipmodel("/home/des/camcat/Data/slipmodels/Tohoku/srcmod/s2011TOHOKU01IDEx", eqfm, nf);
	return 0;
}

int test_nth_index(){

	int Ndim=5;
	int dim[]={2, 4, 3, 5, 5};
	int np=1;
	int *r;

	for (int d=0; d<Ndim; d++) np*=dim[d];
	for (int n=1; n<=np; n++){
		printf("%d:\t", n);
		r=nth_index(n, Ndim, dim);
		for (int d=0; d<Ndim; d++) printf("%d\t",r[d]);
		printf("\n");
	}

	return 1;
}

int test_read_inputfiles(){

	char *file="input/testinput.txt";
	char outname[120], reftime_str[120], fore_template[120], catname[120], focmeccat[120], \
		background_rate_grid[120], slipmodelfile[120], afterslipmodelfile[120];
	struct tm reftime;
	double tstart, tend;
	struct slipmodels_list slip_list;
	int nfm;

	read_inputfile(file, outname, reftime_str, fore_template, catname, focmeccat, background_rate_grid, NULL, NULL,
			slipmodelfile, afterslipmodelfile, NULL, NULL, &reftime, &tstart, &tend, NULL, &nfm);

	printf("outname=%s\n", outname);
	printf("reftime_str=%s\n", reftime_str);
	printf("fore_template=%s\n", fore_template);
	printf("catname=%s\n", catname);
	printf("focmeccat=%s\n", focmeccat);
	printf("background_rate_grid=%s\n", background_rate_grid);
	printf("slipmodelfile=%s\n", slipmodelfile);
	printf("afterslipmodelfile=%s\n", afterslipmodelfile);
	printf("reftime=%d-%d-%dT%d:%d:%dZ\n", reftime.tm_year+1900, reftime.tm_mon+1, reftime.tm_mday, reftime.tm_hour, reftime.tm_min, reftime.tm_sec);
	printf("tstart=%.3lf\n", tstart);
	printf("tend=%.3lf\n", tend);

	read_listslipmodel(slipmodelfile, reftime, &slip_list, 3.0, 0);
	printf("\nslipmodels:\n NSM=%d, is_aft=%d\n", slip_list.NSM, slip_list.is_afterslip);
	for (int n=0; n<slip_list.NSM; n++) printf("t=%.3lf, m=%.3lf, d=%.3lf, file=%s, Nf=%d, no_mod=%d\n", slip_list.tmain[n], slip_list.mmain[n], slip_list.disc[n], slip_list.slipmodels[n], slip_list.Nfaults[n], slip_list.no_slipmodels[n]);

	read_listslipmodel(afterslipmodelfile, reftime, &slip_list, 3.0, 1);
	printf("\nafterslip:\nNSM=%d, is_aft=%d\n", slip_list.NSM, slip_list.is_afterslip);
	for (int n=0; n<slip_list.NSM; n++) printf("t=%.3lf, d=%.3lf, file=%s, Nf=%d, no_mod=%d\n", slip_list.tmain[n], slip_list.disc[0], slip_list.slipmodels[n], slip_list.Nfaults[0], slip_list.no_slipmodels[0]);

	return 0;
}

//int test_background_rate2(){
//
//	char fname[120];
//	FILE *fout;
//	char crust_file[]="input/inCan.dat";
//	char fore_file[]="input/darf_temp.txt";
//	//char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/1863_2011.9.3.dat";
//	//char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/Canterbury_M2_1980_2012.dat";
//	char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/nz_2009_2011_fixed.dat";
//	struct crust crst;
//	struct tm reftime;
//	struct catalog cat;
//	double Mcut=20.0, Mmain=5.0;
//	double dR=50, dZ=50;
//	int ord=2, err;
//	double res=3.0, res_z=1000.0;
//	double t0, dt=365, t00=-20*365;
//	double *weights=NULL;
//
//	sscanf("2010-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
//	reftime.tm_year-=1900;
//	reftime.tm_mon-=1;
//	reftime.tm_isdst=0;
//
//	read_crust(crust_file, fore_file, &crst, res, res_z);
//
//	t0=t00;
//	int c=1;
////	while (t0+dt<=0.0){
////		err=background_rate(cat_file, &crst, reftime, Mcut, Mmain, t0,t0+dt, dR, dZ, ord);
////		if (err) {
////			t0+=dt;
////			continue;
////		}
////		sprintf(fname, "%s/Helm_rate_Darf%d_2nd.dat", testfolder,c);
////		for (int i=1; i<=crst.N_allP; i++) crst.rate0[i]*=crst.r0;
////		print_grid0(fname, crst, crst.rate0);
////
////		sprintf(fname, "%s/Helm_bgrate_Darf%d_2nd.dat", testfolder,c);
////		fout=fopen(fname,"w");
////		fprintf(fout,"Background rate= %.6lf earthquakes/day (M>=%.1lf)\n", crst.r0, crst.mags[1]);
////		fclose(fout);
////
////		t0+=dt;
////		c++;
////	}
//
//	err=background_rate2(cat_file, &crst, reftime, Mcut, Mmain, t00, 0.0, dR, dZ, ord);
//	if (!err){
//		sprintf(fname, "%s/Helm_rate_Darf_dec20092011_2nd.dat", testfolder);
//		for (int i=1; i<=crst.N_allP; i++) crst.rate0[i]*=crst.r0;
//		print_rate(fname, crst, crst.rate0);
//
//		sprintf(fname, "%s/Helm_bgrate_Darf_dec20092011_2nd.dat", testfolder);
//		fout=fopen(fname,"w");
//		fprintf(fout,"Background rate= %.6lf earthquakes/day (M>=%.1lf)\n", crst.r0, crst.mags[1]);
//		fclose(fout);
//	}
//
//	cat.Mc=Mcut;
//	crst.GRmags=dvector(1,1);
//	crst.GRmags[1]=1.0;
//	readZMAP(&cat, NULL, NULL, cat_file, crst, reftime, 0.0, 0.0, -1e30, 0.0, 10, 0.0, 0.0, 0.0, 0.0, 0);
//
//	decluster_catalog(cat, Mmain, &weights, 0);
//
//	sprintf(fname, "%s/Helm_cat_Darf_dec20092011_2nd.dat",testfolder);
//	fout=fopen(fname,"w");
//	for (int i=1; i<=cat.Z; i++) fprintf(fout,"%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",  cat.t[i], cat.lat0[i], cat.lon0[i], cat.depths0[i], cat.mag[i], weights[i]);
//	fclose(fout);
//
//	return 0;
//}

int test_background_rate(){

	char fname[120];
	FILE *fout;
	char crust_file[]="input/inCan.dat";
	char fore_file[]="input/darf_temp.txt";
	char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/nz_2009_2011_fixed.dat";
	struct crust crst;
	struct tm reftime;
	struct catalog cat;
	double Mcut=20, Mmain=7.0;
	double dR=50, dZ=50;
	int ord=1;
	double res=3.0, res_z=1.0;
	double *rates=NULL;

	sscanf("2010-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	reftime.tm_year-=1900;
	reftime.tm_mon-=1;
	reftime.tm_isdst=0;

	//broken (no crust_file anymore)
	//read_crust(crust_file, fore_file, NULL, &crst, res, res_z);
	background_rate(cat_file, &crst, reftime, Mcut, Mmain, -1e30, 1e30,  dR, dZ, ord);

	cat.Mc=Mcut;
	crst.GRmags=dvector(1,1);
	crst.GRmags[1]=1.0;
	readZMAP(&cat, NULL, NULL, cat_file, crst, reftime, 0.0, 0.0, -1e30, 0.0, 10, 0.0, 0.0, 0.0, 0.0, 0);

//	sprintf(fname, "%s/Helm_cat_Darf_past_3d_MC3.0.dat",testfolder);
//	fout=fopen(fname,"w");
//	for (int i=1; i<=cat.Z; i++) fprintf(fout,"%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", cat.lat0[i], cat.lon0[i], cat.depths0[i], cat.t[i], cat.mag[i]);
//	fclose(fout);

	sprintf(fname, "%s/Background_Darf_hr.dat", testfolder);
	for (int i=1; i<=crst.N_allP; i++) crst.rate0[i]*=crst.r0;
	print_rate(fname, crst, 3.0, crst.rate0);

	convert_geometry(crst, crst.rate0, &rates, 1, 0);
	sprintf(fname, "%s/Background_Darf.dat", testfolder);
	csep_forecast(fname, crst, rates, 0);

//	sprintf(fname, "%s/Helm_rate_Darf_3d_MC3.0.dat", testfolder);
//	fout=fopen(fname,"w");
//	for (int i=1; i<=crst.N_allP; i++) {
//		fprintf(fout,"%.5lf\t%.5lf\t%.5lf\t%.5e\n",crst.lat[i], crst.lon[i], crst.depth[i], crst.rate0[i]);
//	}
//	fclose(fout);
//
//	sprintf(fname, "%s/Helm_bgrate_Darf_3d_MC3.0.dat", testfolder);
//	fout=fopen(fname,"w");
//	fprintf(fout,"Background rate= %.6lf earthquakes/day (M>=%.1lf)\n", crst.r0, crst.mags[1]);
//	fclose(fout);

	return 0;
}

int test_decluster_catalog(){

	char fname[120];
	double Mmain=6.0;
	double *time_missing=NULL;
	FILE *fout;
	char crust_file[]="input/Tohoku_simple_vert.inp";
	char fore_file[]="input/tohoku_template.dat";
	char cat_file[]="/home/des/camcat/Data/Catalogs/Others/jma_cat_2010_2013_update20130329_sel_2.5.dat";
	struct catalog cat;
	struct crust crst;
	struct tm tt;
	time_t t;
	int *dec;
	double *weights;
	double *rate, r;
	double res=6.0;

	time(&t);
	tt=*(localtime(&t));

	//broken (no crust_file anymore)
	//read_crust(crust_file, fore_file, NULL, &crst, res, 100.0);
	cat.Mc=0.0;
	readZMAP(&cat, NULL, NULL, cat_file, crst, tt, 0.0, 0.0, -1e30, 1e30, 10, 0.0, 0.0, 0.0, 0.0, 0);
	weights=dvector(1,cat.Z);


	//old declustering method (rescales rates after calculating them):
	dec=decluster_catalog_rescalegrid(cat, crst, Mmain, &time_missing, 0);

	for (int i=1; i<=cat.Z; i++) weights[i]=(double) dec[i];

	rate=Helmstetter_cat(cat, crst, weights, 2);

	sprintf(fname, "%s/Helm_rate_dec6.0_2nd.dat", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=crst.nLat*crst.nLon; i++) fprintf(fout,"%.5lf\t%.5lf\t%.5lf\n",crst.lat[i], crst.lon[i],rate[i]);
	fclose(fout);

	sprintf(fname, "%s/Helm_rate_dec_rescaled6.0_2nd.dat", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=crst.nLat*crst.nLon; i++) {
		r=(cat.tend-cat.tstart)/(cat.tend-cat.tstart-time_missing[i]);
		fprintf(fout,"%.5lf\t%.5lf\t%.5lf\n",crst.lat[i], crst.lon[i],r*rate[i]);
	}
	fclose(fout);

	//new declustering method (weights earthquakes):
	decluster_catalog(cat, Mmain, &weights, 0);

	sprintf(fname, "%s/Dec_cat6.0.dat",testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=cat.Z; i++) fprintf(fout,"%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", cat.lat0[i], cat.lon0[i], cat.depths0[i], cat.mag[i], cat.t[i], weights[i]);
	fclose(fout);

	rate=Helmstetter_cat(cat, crst, weights, 2);
	sprintf(fname, "%s/Helm_rate_dec_rescaled_new6.0_2nd.dat", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=crst.nLat*crst.nLon; i++) fprintf(fout,"%.5lf\t%.5lf\t%.5lf\n",crst.lat[i], crst.lon[i],rate[i]);
	fclose(fout);

	return 0;

}

int test_fit_depth(){

	char fname[120];
	FILE *fout;
	char crust_file[]="input/Tohoku_simple_vert.inp";
	char fore_file[]="input/tohoku_template.dat";
	char cat_file[]="/home/des/camcat/Data/Catalogs/Others/jma_cat_2010_2013_update20130329_sel_2.5.dat";
	struct catalog cat;
	struct crust crst;
	struct tm tt;
	time_t t;
	double *rate;
	double res=6.0;
	double *depths, *p, *verr;

	time(&t);
	tt=*(localtime(&t));

	//broken (no crust_file anymore)
	//read_crust(crust_file, fore_file,NULL,  &crst, 100.0, 1.0);
//	gridPMax=crst.N_allP;
	gridPMax=1000;
	cat.Mc=0.0;
	readZMAP(&cat, NULL, NULL, cat_file, crst, tt, 0.0, 0.0, -1e30, 1e30, 10, 0.0, 0.0, 0.0, 0.0, 0);

	depths=dvector(1,crst.nD);
	for (int i=1; i<=crst.nD; i++) {
		depths[i]=crst.depth[1+(i-1)*crst.nLat*crst.nLon];
	}

	verr=dvector(1,cat.Z);
	for (int i=1; i<=cat.Z; i++) verr[i]=10;

	p=fit_depth(depths, crst.ddepth, crst.nD, cat.depths0, verr, NULL, cat.Z);

	sprintf(fname, "%s/Helm_catV3a.dat",testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=cat.Z; i++) fprintf(fout,"%.3lf\n", cat.depths0[i]);
	fclose(fout);

	sprintf(fname, "%s/Helm_rateV3a.dat", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=crst.nD; i++) fprintf(fout,"%.5lf\t%.5lf\n",depths[i], p[i]);
	fclose(fout);

	return(0);
}

int test_Helmstetter_cat(){

	char fname[120];
	FILE *fout;
	//Tohoku:
//	char crust_file[]="input/Tohoku_simple_vert.inp";
//	char fore_file[]="input/tohoku_template.dat";
//	char cat_file[]="/home/des/camcat/Data/Catalogs/Others/jma_cat_2010_2013_update20130329_sel_2.5.dat";
	//Darfield:
	char crust_file[]="input/inCan.dat";
	char fore_file[]="input/darf_temp.txt";
	//char cat_file[]="input/quake_recent.dat";
	char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new zeland/nz_2009_2011_fixed.dat";
	struct catalog cat;
	struct crust crst;
	struct tm reftime;
	double *rate;
	double res=3.0;

	sscanf("2010-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	reftime.tm_year-=1900;
	reftime.tm_mon-=1;
	reftime.tm_isdst=0;

	//broken (no crust_file anymore)
	//read_crust(crust_file, fore_file, NULL, &crst, res, 100.0);
//	gridPMax=crst.N_allP;
	gridPMax=1000;
	cat.Mc=0.0;
	readZMAP(&cat, NULL, NULL, cat_file, crst, reftime, 0.0, 0.0, -1e30, 0.0, 10, 0.0, 0.0, 0.0, 0.0, 0);

	rate=Helmstetter_cat(cat, crst, NULL, 2);

	sprintf(fname, "%s/Helm_cat3_Darf_past_2nd.dat",testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=cat.Z; i++) fprintf(fout,"%.3lf\t%.3lf\t%.3lf\t%.3lf\n", cat.lat0[i], cat.lon0[i], cat.t[i], cat.mag[i]);
	fclose(fout);

	sprintf(fname, "%s/Helm_rate3_Darf_past_2nd.dat", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=crst.nLat*crst.nLon; i++) fprintf(fout,"%.5lf\t%.5lf\t%.5lf\n",crst.lat[i], crst.lon[i],rate[i]);
	fclose(fout);

	return 0;

}

int test_Helmstetter(){

	char fname[120];
	FILE *fout;
	long seed=-19329935;
	int ne=100, N=100, N3=N*N;
	double *xe, *ye, *err;
	double *xs, *ys, *rate;
	double dx=1.0/(1.0*N), dy=dx;

	xs=dvector(1,N3);
	ys=dvector(1,N3);
	err=dvector(1,N3);
	for (int i=1; i<=N; i++) {
		for (int j=1; j<=N; j++) {
			xs[(i-1)*N+j]=(i-1.0)*dx;
			ys[(i-1)*N+j]=(j-1.0)*dy;
			err[(i-1)*N+j]=1000.0;
		}
	}

	xe=dvector(1,N);
	ye=dvector(1,N);
	for (int i=1; i<=N; i++) {
		xe[i]=ran1(&seed);
		seed*=1.0;
		ye[i]=ran1(&seed);
	}

	rate=Helmstetter(xs, ys, dx, dy, N3, xe, ye, err, NULL, ne, 1);

	sprintf(fname, "%s/Helm_cat.dat",testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=N; i++) fprintf(fout,"%.3lf\t%.3lf\n", xe[i], ye[i]);
	fclose(fout);

	sprintf(fname, "%s/Helm_rate.dat", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=N3; i++) fprintf(fout,"%.5lf\t%.5lf\t%.5lf\n",xs[i], ys[i],rate[i]);
	fclose(fout);

	return 0;
}

int test_mysort(){

	int N=10;
	double xs[]={5.0,6.0,3.0,4.0,9.0,10.0,0.0,1.0,2.0,6.0,8.0};
	int *ind=NULL;
	double *arr=NULL;

	mysort(N, xs, &ind, &arr);

	for (int i=1; i<=N; i++) printf("%.0lf\t", xs[ind[i]]);

	return 0;
}

int test_all_2ndnearestneighbours(){

	char fname[120];
	FILE *fout;
	long seed=-19329935;
	int *pts=NULL;
	double *xs, *ys, *dist=NULL;
	int N=100, n;

	sprintf(fname, "%s/nearest_neighbour_2nd_sparse.dat",testfolder);

	xs=dvector(1,N);
	ys=dvector(1,N);

	for (int i=1; i<=N; i++) {
		xs[i]=ran1(&seed);
		seed*=1.0;
		ys[i]=ran1(&seed);
	}

	n=all_2ndnearestneighbours(xs, ys, N, &pts, &dist);
	printf("N=%d, n=%d\n", N, n);

	fout=fopen(fname,"w");
	for (int i=1; i<=N; i++) fprintf(fout,"%.3lf\t%.3lf\t%d\t%.5e\n", xs[i], ys[i], pts[i], dist[i]);
	fclose(fout);

	return(0);

}

int test_all_nearestneighbours(){
// test gave 221 operations with N=100, <28000 for N=10000.

	char fname[120];
	FILE *fout;
	long seed=-19329935;
	int *pts=NULL;
	double *xs, *ys, *dist=NULL;
	int N=100, n;

	sprintf(fname, "%s/nearest_neighbour_sparse.dat",testfolder);

	xs=dvector(1,N);
	ys=dvector(1,N);

	for (int i=1; i<=N; i++) {
		xs[i]=ran1(&seed);
		seed*=1.0;
		ys[i]=ran1(&seed);
	}

	n=all_nearestneighbours(xs, ys, N, &pts, &dist);
	printf("N=%d, n=%d\n", N, n);

	fout=fopen(fname,"w");
	for (int i=1; i<=N; i++) fprintf(fout,"%.3lf\t%.3lf\t%d\t%.5e\n", xs[i], ys[i], pts[i], dist[i]);
	fclose(fout);

	return(0);

}

int test_find_gridpoints_exact(){
	/*
	 * output produced by this function (so far):
	 * find_gridpointsXX3: uses SD=0.2 (set inside functions find_gidpointXX). xe=ye=0.5
	 * find_gridpointsXX4: uses sd given below. xe=ye=0.55
	 */

	char fname[120];
	FILE *fout, *fout0;
	int N=8, N3=N*N;
	int np, *points;
	double *weights, *weights0;
	double *xs, *ys, *zs, *dAs;
	double dx=1.0/(1.0*N), dy=dx;
	double sd=0.001;
	double xe=0.55, ye=0.55, ze=0;

	xs=dvector(1,N3);
	ys=dvector(1,N3);
	zs=dvector(1,N3);
	dAs=dvector(1,N3);
	points=ivector(1,N3);
	weights=dvector(1,N3);
	weights0=dvector(1,N3);

	for (int i=1; i<=N; i++) {
		for (int j=1; j<=N; j++) {
			xs[(i-1)*N+j]=(i-1.0)*dx;
			ys[(i-1)*N+j]=(j-1.0)*dy;
			zs[(i-1)*N+j]=0.0;
			dAs[(i-1)*N+j]=1.0;
		}
	}

	//find_gridpoints(ys, xs, dAs, zs, N3, N3, ye, xe, sd, ze, sd, 10000, &np, points, weights0, 0);
	find_gridpoints_exact(ys, xs, zs, dx, dy, 0.0, N3, N3, ye, xe, sd, ze, sd, 10000, &np, points, weights, 0, 0);

	sprintf(fname, "%s/find_gridpointsapprox4",testfolder);
	fout0=fopen(fname,"w");
	sprintf(fname, "%s/find_gridpointsexact4", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=N; i++){
		for (int j=1; j<=N; j++){
			fprintf(fout,"%.5lf\t",weights[(i-1)*N+j]);
			fprintf(fout0,"%.5lf\t",weights0[(i-1)*N+j]);
		}
		fprintf(fout,"\n");
		fprintf(fout0,"\n");
	}
	fclose(fout);
	fclose(fout0);

	free_dvector(xs, 1, N3);
	free_dvector(ys, 1, N3);
	free_dvector(zs, 1, N3);
	free_dvector(dAs, 1, N3);
	free_dvector(weights0, 1, N3);
	free_dvector(weights, 1, N3);
	free_ivector(points, 1, N3);

	return 0;

}

int test_exact_prob(){

	char fname[120];
	FILE *fout, *fout0;
	int N=8;
	double *xs, *ys;
	double dx=1.0/(1.0*N), dy=dx;
	double sd=0.05;
	double xe=0.5, ye=0.5;
	double **p, p_tot=0;
	double **p0, p0_tot=0;

	xs=dvector(1,N);
	ys=xs;
	p=dmatrix(1,N,1,N);
	p0=dmatrix(1,N,1,N);

	for (int i=1; i<=N; i++) xs[i]=(i-1.0)*dx;

	for (int i=1; i<=N; i++){
		for (int j=1; j<=N; j++){
			p0[i][j]=exp(-pow(xs[i]-xe,2)/(2*pow(sd,2)))*exp(-pow(ys[j]-ye,2)/(2*pow(sd,2)));
			p0_tot+=p0[i][j];
			p[i][j]=exact_prob(xs[i]-xe, ys[j]-ye, 0.0, dx, dy, 0.0, sd, sd, 0.0, 0);
			p_tot+=p[i][j];
		}
	}

	sprintf(fname, "%s/find_gridpointsapprox2",testfolder);
	fout0=fopen(fname,"w");
	sprintf(fname, "%s/find_gridpointsexact2", testfolder);
	fout=fopen(fname,"w");
	for (int i=1; i<=N; i++){
		for (int j=1; j<=N; j++){
			fprintf(fout,"%.5lf\t",p[i][j]/p_tot);
			fprintf(fout0,"%.5lf\t",p0[i][j]/p0_tot);
		}
		fprintf(fout,"\n");
		fprintf(fout0,"\n");
	}
	fclose(fout);
	fclose(fout0);

	free_dvector(xs, 1, N);
	free_dmatrix(p, 1, N,1,N);
	free_dmatrix(p0, 1, N,1,N);

	return 0;

}

int test_readZMAP_tw(){

	struct catalog cat;
	struct crust crst;
	struct tm reftime;
	int NT;
	char fname[120];
	char *file="/home/des/camcat/Data/Catalogs/Others/jma_cat_2010_2013_update20130329_sel_2.5.dat";
	char *crust_file="/home/des/camcat/Code/CRS_1.0/INPUT/Tohoku_new.inp";
	char *fore_template="input/tohoku_template_sparse.dat";
	double t0=0.0, t1=300, tw=1.0;
	double Mmain=6.8;
	double res=100;

	gridPMax=1000;
//	2011-03-11T14:46:18Z

	setenv("TZ", "UTC", 1);
	reftime.tm_year=111;
	reftime.tm_mon=2;
	reftime.tm_isdst=0;
	reftime.tm_mday=03;
	reftime.tm_hour=14;
	reftime.tm_min=46;
	reftime.tm_sec=18;

	//broken (no crust_file anymore)
	//read_crust(crust_file, fore_template, NULL,  &crst, res, res);
	readZMAP (&cat, NULL, &NT, file, crst, reftime, t0,t1, t0, t1, Mmain, 0.0, 0, 0, 1e5, 1);
	sprintf(fname, "%s/cat_notw.dat", testfolder);
	print_cat(fname, cat);
	readZMAP (&cat, NULL, &NT, file, crst, reftime, t0,t1, t0, t1, Mmain, tw, 0, 0, 1e5, 1);
	sprintf(fname, "%s/cat_tw.dat", testfolder);
	print_cat(fname, cat);

	return 0;
}

int test_forecast_stepG2_new(){
//tests forecast after stress steps (and potentially, also afterslip).

	char name[]="1_2aslipeqks";	//file name.
	int NP=10, Neq=3, NTS=2, n_samples=3000;
	struct catalog cat;
	struct pscmp *DCFS;
	double times[NTS-1];
	double **cmpdata;	//afterslip.
	double Asig=1000, ta=1000;
	double cmb_step=20000;
	double *gamma0;
	double *Nend=dvector(1,n_samples);
	double *Rend=dvector(1,n_samples);
	double t0=0, t1=300;
	double dt=(t1-t0)/n_samples;
	int *points;
	char fname[120];
	double f=0.5;
	FILE *fout;

	points=ivector(1,NP);
	for (int j=1; j<=NP; j++) points[j]=j;

	for (int i=0; i<NTS; i++) times[i]= t0-0.001+i*(t1+0.001-t0)/(NTS-1);

	cmpdata=dmatrix(0,NTS,1,NP);
	for (int i=0; i<=NTS; i++){
		for (int j=1; j<=NP; j++){
			cmpdata[i][j]=f*cmb_step*(log(2+i)-log(1+i))/log(2+NTS);
		}
	}

	DCFS=pscmp_array(0,Neq);
	for (int i=0; i<Neq; i++){
		DCFS[i].t=t0+i*(t1-t0)/10;
		DCFS[i].m=3.0;
		DCFS[i].nsel=NP;
		DCFS[i].which_pts=points;
		DCFS[i].cmb=dvector(1,NP);
		for (int j=1; j<=NP; j++) DCFS[i].cmb[j]=cmb_step;
	}

	cat.Z=0;

	gamma0=dvector(1,NP);
	for (int j=1; j<=NP; j++) gamma0[j]=Asig/ta;

	//for (int t=1; t<=n_samples; t++) forecast_stepG2_new(cat, times, cmpdata, DCFS, t0+dt*(t-1), t0+dt*t, Asig, ta, points, NULL, Nend+t, Rend+t, NP, NTS, Neq, gamma0, NULL, 1);

	sprintf(fname, "%s/forecast_stepG2_new3/%s.readme", testfolder, name);
	fout=fopen(fname,"w");
	fprintf(fout, "Asig\t\tta\tcmb_step\t\tafterslip_cmb_step\tt0\tt1\tdt\n");
	fprintf(fout, "%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.5lf\n", Asig, ta, cmb_step, f*cmb_step, t0, t1, fmin((t1+0.001-t0)/(NTS-1), (t1-t0)/n_samples));
	fprintf(fout, "Output format: \ntime(days)	cumu no.\trate(1/days)\n");
	fclose(fout);

	sprintf(fname, "%s/forecast_stepG2_new3/%s.dat", testfolder, name);
	fout=fopen(fname,"w");
	for (int i=1; i<=n_samples; i++) fprintf(fout, "%lf\t%lf\t%lf\n", t1+dt*i, Nend[i], Rend[i]);
	fclose(fout);

	printf("Done!\n");
	return 0;

}

int test_readtxttemplate(){

	int nmag, nlat, nlon, ndep, ng;
	double dlat, dlon, ddep;
	double *lats, *lons, *deps;
	double minlat, maxlat, minlon, maxlon, mindep, maxdep, m0, m1;

	read_csep_template("input/darf_temp.txt", &nmag, &nlat, &nlon, &ndep, &ng, &dlat, &dlon, &ddep, NULL, &lats, &lons, &deps, 0, &minlat, &maxlat, &minlon, &maxlon, &mindep, &maxdep, &m0, &m1, NULL);

	printf("ng=%d, nlat=%d, nlon=%d, ndep=%d, nmag=%d\n", ng, nlat, nlon, ndep, nmag);
	printf("dlat=%lf, dlon=%lf, ddep=%lf\n", dlat, dlon, ddep);
	printf("minlat=%lf, minlon=%lf, mindep=%lf, minmag=%lf\n", minlat, minlon, mindep, m0);
	printf("maxlat=%lf, maxlon=%lf, maxdep=%lf, maxmag=%lf\n", maxlat, maxlon, maxdep, m1);

	return 0;
}

int test_pointers(){

	double *a, *b;
	a=dvector(1,10);

	for (int i=1; i<=10; i++) a[i]=i;
	b=a;

	a=dvector(1,10);

	for (int i=1; i<=10; i++) a[i]=i*i;

	printf("a=[");
	for (int i=1; i<=10; i++) printf("%lf, ", a[i]);
	printf("]\n b=[");
	for (int i=1; i<=10; i++) printf("%lf, ", b[i]);
	printf("]\n");
}

int fun_with_static(){

	static double c=0;
	c+=1;
	printf("c=%lf\n",c);
	sleep(1);
	return 0;
}

int test_staticparallel(){

	#pragma omp parallel for
	for (int i=1; i<=100; i++) {
		printf("%d\t", omp_get_thread_num());
		fun_with_static();
	}
	return 0;
}

int test_allOkada_simple_multiplerec(){
/* calculates a DCFS field from eqkfm structure
 * to visualize in Paraview, use following commands:
 * for i in $(seq 0 9); do /home/des/camcat/Code/Scripts/slipmodel_sep.sh slipmodel$i.dat;
 * python /home/des/camcat/Code/Scripts/Slipmodel2vtk.py slipmodel$(echo $i)_ Tohoku;
 * python /home/des/camcat/Code/Scripts/Foremap2vtk.py okadaDCFS$i.dat Tohoku; done
 */

	int N=6;
	long seed=-19329935;
	double res=6.0;
	struct eqkfm eqfm;
	struct crust crst;
	FILE *fout;
	char fname[120], fname0[120];
	struct pscmp dcfs;
	struct Coeff_LinkList coeffs;
	float **resCoeff_st, **resCoeff_di;
	double lat0=-38.5, lon0=142.0, dep0=0.0;
	int Nlat=20, Nlon=40, Ndep=10, NP=Nlat*Nlon*Ndep;
	double *latgrid, *longrid, *depgrid;
	int *allp, nofm=1;
	double dlat=1.2, dlon=1.2, ddep=4.0;
	double strikes[nofm], dips[nofm], rakes[nofm];
	time_t time1, time2;

	latgrid=dvector(1,NP);
	longrid=dvector(1,NP);
	depgrid=dvector(1,NP);
	allp=ivector(1,NP);
	crst.x=dvector(1,NP);
	crst.y=dvector(1,NP);

	double strikes0[]={90, 90, 90, 90, 90, 90};
	double dips0[]={90, 90, 45, 15, 15, 15};
	double rakes0[]={0, 180, 45, 90, 270, 45};

	for (int i=1; i<=NP; i++) allp[i]=i;

	for (int i=1; i<=Nlat; i++){
		for (int j=1; j<=Nlon; j++){
			for (int k=1; k<=Ndep; k++){
				latgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lat0-0.5*dlat+(dlat/Nlat)*(i-0.5);
				longrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lon0-0.5*dlon+(dlon/Nlon)*(j-0.5);
				depgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= (ddep/Ndep)*(k-0.5);
			}
		}
	}

	//broken: function does not exist anymore
	//read_pscmp_crust("/home/des/camcat/Code/CRS_2.01/input/Tohoku_simple_new.inp",&crst);
	sprintf(fname0, "%s/okada3/fm_slipmodels_summary.txt",testfolder);
	fout=fopen(fname0,"w");

	fprintf(fout, "strike\t dip \t rake \t magnitude \t resolution \n");
	eqfm.lat=lat0;
	eqfm.lon=lon0;
	eqfm.whichfm=1;
	eqfm.taper=ivector(1,4);
	for (int i=1; i<=3; i++) eqfm.taper[i]=1;
	eqfm.taper[4]=0;

	dcfs.NF=1;
	crst.fmzone=ivector(1,NP);
	for (int i=1; i<=NP; i++) crst.fmzone[i]=(int) (nofm-1)*ran1(&seed);
	crst.nLat=Nlat;
	crst.nLon=Nlon;
	crst.nD=Ndep;
	crst.lat0=lat0;
	crst.lon0=lon0;
	crst.lat=latgrid;
	crst.lon=longrid;
	crst.depth=depgrid;
	crst.dlat=dlat/Nlat;
	crst.dlon=dlon/Nlon;
	crst.ddepth=ddep/Ndep;
	crst.nofmzones=nofm;
	for (int i=1; i<=NP; i++) latlon2localcartesian(latgrid[i], longrid[i], lat0, lon0, crst.y+i, crst.x+i);
	latlon2localcartesian(eqfm.lat, eqfm.lon, crst.lat0, crst.lon0, &(eqfm.y), &(eqfm.x));


	time(&time1);

	for (int n=0; n<N; n++){

		eqfm.depth=dep0;
		eqfm.str1=strikes0[n];
		eqfm.dip1=dips0[n];
		eqfm.rake1=rakes0[n];
		eqfm.whichfm=1;
		eqfm.mag=6.3;
		//res=0.5+ran1(&seed);

		for (int fm=0; fm<=nofm; fm++) {
			strikes[fm]=eqfm.str1;
			dips[fm]=eqfm.dip1;
			rakes[fm]=eqfm.rake1;
		}

		fprintf(fout, "%.0lf\t%.0lf\t%.0lf\t%.1lf\t%.1lf\n", eqfm.str1, eqfm.dip1, eqfm.rake1, eqfm.mag, res);
		focmec2slipmodel(crst, &eqfm, res, 0, 0);
		dcfs.nsel=eqfm.nsel=NP;
		dcfs.which_pts=eqfm.selpoints=allp;
		printf("%d\n",dcfs.nsel);
		dcfs.S=d3tensor(1,dcfs.nsel, 1,3,1,3);
		dcfs.cmb=dvector(1,dcfs.nsel);

		// 1) using okadaDCFS:

		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;
		okadaDCFS(dcfs, &eqfm, 1, crst, strikes, dips, 1);
		resolve_DCFS(dcfs, crst, strikes, dips, rakes, 0);
		sprintf(fname, "%s/okada3/okadaDCFSfixA%d.dat",testfolder, n);	//all wrong
		//sprintf(fname, "%s/okadaCoeff4%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaDCFS(dcfs, &eqfm, 1, crst, NULL, NULL, 1);
		resolve_DCFS(dcfs, crst, strikes, dips, NULL, 1);
		sprintf(fname, "%s/okada3/okadaDCFSfix_optrakeB%d.dat",testfolder, n);	//0,2,3 wrong
		//sprintf(fname, "%s/okadaCoeff5%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		resolve_DCFS(dcfs, crst, strikes, dips, rakes, 0);
		sprintf(fname, "%s/okada3/okadaDCFSfixB%d.dat",testfolder, n);	//all wrong
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;
		DCFScmbopt(&dcfs, 0, crst);
		sprintf(fname, "%s/okada3/okadaDCFSopt%d.dat",testfolder, n);	//all right
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

//		//2) using okadaCoeff:
//
		coeffs.NgridT=dcfs.nsel;
		coeffs.NP=eqfm.np_di*eqfm.np_st;
		okadaCoeff(&(coeffs.Coeffs_st), &(coeffs.Coeffs_dip), &eqfm, 1, crst, latgrid, longrid, depgrid);

		okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst,  strikes, dips, 0);
		sprintf(fname, "%s/okada3/CRS_2.01%d.dat",testfolder, n);	//0,2,3 wrong
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

	//	okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst, 0.0, 0.0, 1);
		resolve_DCFS(dcfs, crst, strikes, dips, NULL, 1);
		sprintf(fname, "%s/okada3/okadaCoeff_fix_optrakeB%d.dat",testfolder, n);	//0,2,3 wrong
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst, NULL, NULL, 1);
		resolve_DCFS(dcfs, crst, strikes, dips, rakes, 0);
		sprintf(fname, "%s/okada3/okadaCoeff_fix_A%d.dat",testfolder, n);	//all wrong
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaCoeff_resolve(coeffs, &resCoeff_st, &resCoeff_di, crst, strikes, dips, rakes);
		resolvedCoeff2DCFS(resCoeff_st, resCoeff_di, dcfs, &eqfm, crst);
		sprintf(fname, "%s/okada3/okadaCoeff_fix_B%d.dat",testfolder, n);	//all wrong
		//sprintf(fname, "%s/okadaCoeff4%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);

		free_dvector(dcfs.cmb,1,dcfs.nsel);
		free_d3tensor(dcfs.S,1,dcfs.nsel,1,3,1,3);
		free_f3tensor(coeffs.Coeffs_dip,1,coeffs.NP, 1, coeffs.NgridT, 1,6);
		free_f3tensor(coeffs.Coeffs_st,1,coeffs.NP, 1, coeffs.NgridT, 1,6);
		free_matrix(resCoeff_st, 1,coeffs.NP, 1, coeffs.NgridT);
		free_matrix(resCoeff_di, 1,coeffs.NP, 1, coeffs.NgridT);
	}
	fclose(fout);

	time(&time2);
	printf("Execution time: %.3lf seconds.\n", difftime(time2, time1));

	free_dvector(latgrid,1,NP);
	free_dvector(longrid,1,NP);
	free_dvector(depgrid,1,NP);
	free_dvector(crst.x,1,NP);
	free_dvector(crst.y,1,NP);
	free_ivector(allp,1,NP);

	return 0;

}

int test_allOkada_simple(){
/* calculates a DCFS field from eqkfm structure
 * to visualize in Paraview, use following commands:
 * for i in $(seq 0 9); do /home/des/camcat/Code/Scripts/slipmodel_sep.sh slipmodel$i.dat;
 * python /home/des/camcat/Code/Scripts/Slipmodel2vtk.py slipmodel$(echo $i)_ Tohoku;
 * python /home/des/camcat/Code/Scripts/Foremap2vtk.py okadaDCFS$i.dat Tohoku; done
 */

	int N=6;
	//long seed=-19329935;
	double res=6.0;
	struct eqkfm eqfm;
	struct crust crst;
	FILE *fout;
	char fname[120], fname0[120];
	struct pscmp dcfs;
	struct Coeff_LinkList coeffs;
	float **resCoeff_st, **resCoeff_di;
	double lat0=-38.5, lon0=142.0, dep0=0.0;
	int Nlat=20, Nlon=40, Ndep=10, NP=Nlat*Nlon*Ndep;
	double *latgrid, *longrid, *depgrid;
	int *allp;
	double dlat=1.2, dlon=1.2, ddep=4.0;
	time_t time1, time2;

	latgrid=dvector(1,NP);
	longrid=dvector(1,NP);
	depgrid=dvector(1,NP);
	allp=ivector(1,NP);
	crst.x=dvector(1,NP);
	crst.y=dvector(1,NP);

	double strikes[]={90, 90, 90, 90, 90, 90};
	double dips[]={90, 90, 45, 15, 15, 15};
	double rakes[]={0, 180, 45, 90, 270, 45};

	for (int i=1; i<=NP; i++) allp[i]=i;

	for (int i=1; i<=Nlat; i++){
		for (int j=1; j<=Nlon; j++){
			for (int k=1; k<=Ndep; k++){
				latgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lat0-0.5*dlat+(dlat/Nlat)*(i-0.5);
				longrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lon0-0.5*dlon+(dlon/Nlon)*(j-0.5);
				depgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= (ddep/Ndep)*(k-0.5);
			}
		}
	}

	//broken: function does not exist anymore
	//read_pscmp_crust("/home/des/camcat/Code/CRS_2.01/input/Tohoku_simple_new.inp",&crst);
	sprintf(fname0, "%s/okada2/fm_slipmodels_summary.txt",testfolder);
	fout=fopen(fname0,"w");

	fprintf(fout, "strike\t dip \t rake \t magnitude \t resolution \n");
	eqfm.lat=lat0;
	eqfm.lon=lon0;
	eqfm.whichfm=1;
	eqfm.taper=ivector(1,4);
	for (int i=1; i<=3; i++) eqfm.taper[i]=1;
	eqfm.taper[4]=0;

	dcfs.NF=1;
	crst.nLat=Nlat;
	crst.nLon=Nlon;
	crst.nD=Ndep;
	crst.lat0=lat0;
	crst.lon0=lon0;
	crst.lat=latgrid;
	crst.lon=longrid;
	crst.depth=depgrid;
	crst.dlat=dlat/Nlat;
	crst.dlon=dlon/Nlon;
	crst.ddepth=ddep/Ndep;
	for (int i=1; i<=NP; i++) latlon2localcartesian(latgrid[i], longrid[i], lat0, lon0, crst.y+i, crst.x+i);
	latlon2localcartesian(eqfm.lat, eqfm.lon, crst.lat0, crst.lon0, &(eqfm.y), &(eqfm.x));


	time(&time1);

	for (int n=0; n<N; n++){
		eqfm.depth=dep0;
		eqfm.str1=strikes[n];
		eqfm.dip1=dips[n];
		eqfm.rake1=rakes[n];
		eqfm.whichfm=1;
		eqfm.mag=6.3;
		//res=0.5+ran1(&seed);

		fprintf(fout, "%.0lf\t%.0lf\t%.0lf\t%.1lf\t%.1lf\n", eqfm.str1, eqfm.dip1, eqfm.rake1, eqfm.mag, res);
		focmec2slipmodel(crst, &eqfm, res, 0, 0);
		dcfs.nsel=eqfm.nsel=NP;
		dcfs.which_pts=eqfm.selpoints=allp;
		printf("%d\n",dcfs.nsel);
		dcfs.S=d3tensor(1,dcfs.nsel, 1,3,1,3);
		dcfs.cmb=dvector(1,dcfs.nsel);

		// 1) using okadaDCFS:

		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;
		okadaDCFS(dcfs, &eqfm, 1, crst, &eqfm.str1, &eqfm.dip1, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, &(eqfm.rake1), 0);
		sprintf(fname, "%s/okada2/okadaDCFSfixA%d.dat",testfolder, n);
		//sprintf(fname, "%s/okadaCoeff4%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaDCFS(dcfs, &eqfm, 1, crst, NULL, NULL, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, NULL, 1);
		sprintf(fname, "%s/okada2/okadaDCFSfix_optrakeB%d.dat",testfolder, n);
		//sprintf(fname, "%s/okadaCoeff5%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, &(eqfm.rake1), 0);
		sprintf(fname, "%s/okada2/okadaDCFSfixB%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;
		DCFScmbopt(&dcfs, 0, crst);
		sprintf(fname, "%s/okada2/okadaDCFSopt%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

//		//2) using okadaCoeff:
//
		coeffs.NgridT=dcfs.nsel;
		coeffs.NP=eqfm.np_di*eqfm.np_st;
		okadaCoeff(&(coeffs.Coeffs_st), &(coeffs.Coeffs_dip), &eqfm, 1, crst, latgrid, longrid, depgrid);

		okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst,  &eqfm.str1, &eqfm.dip1, 0);
		sprintf(fname, "%s/okada2/CRS_2.01%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

	//	okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst, 0.0, 0.0, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, NULL, 1);
		sprintf(fname, "%s/okada2/okadaCoeff_fix_optrakeB%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst, NULL, NULL, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, &(eqfm.rake1),  0);
		sprintf(fname, "%s/okada2/okadaCoeff_fix_A%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaCoeff_resolve(coeffs, &resCoeff_st, &resCoeff_di, crst, &eqfm.str1, &eqfm.dip1, &eqfm.rake1);
		resolvedCoeff2DCFS(resCoeff_st, resCoeff_di, dcfs, &eqfm, crst);
		sprintf(fname, "%s/okada2/okadaCoeff_fix_B%d.dat",testfolder, n);
		//sprintf(fname, "%s/okadaCoeff4%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);

		free_dvector(dcfs.cmb,1,dcfs.nsel);
		free_d3tensor(dcfs.S,1,dcfs.nsel,1,3,1,3);
		free_f3tensor(coeffs.Coeffs_dip,1,coeffs.NP, 1, coeffs.NgridT, 1,6);
		free_f3tensor(coeffs.Coeffs_st,1,coeffs.NP, 1, coeffs.NgridT, 1,6);
		free_matrix(resCoeff_st, 1,coeffs.NP, 1, coeffs.NgridT);
		free_matrix(resCoeff_di, 1,coeffs.NP, 1, coeffs.NgridT);
	}
	fclose(fout);

	time(&time2);
	printf("Execution time: %.3lf seconds.\n", difftime(time2, time1));

	free_dvector(latgrid,1,NP);
	free_dvector(longrid,1,NP);
	free_dvector(depgrid,1,NP);
	free_dvector(crst.x,1,NP);
	free_dvector(crst.y,1,NP);
	free_ivector(allp,1,NP);

	return 0;

}

int test_allOkada(){
/* test all the ways to calculate DCFS field from eqkfm structure (functions in okadaDCFS.c)
 * to visualize in Paraview, use following commands:
 * for i in $(seq 0 9); do /home/des/camcat/Code/Scripts/slipmodel_sep.sh slipmodel$i.dat;
 * python /home/des/camcat/Code/Scripts/Slipmodel2vtk.py slipmodel$(echo $i)_ Darfield;
 * python /home/des/camcat/Code/Scripts/Foremap2vtk.py okadaDCFS$i.dat Darfield; done
 */

	int N=3;
	//long seed=-19329935;
	double res=6.0;
	struct eqkfm eqfm;
	struct crust crst;
	FILE *fout;
	char fname[120], fname0[120];
	struct pscmp dcfs;
	struct Coeff_LinkList coeffs;
	float **resCoeff_st, **resCoeff_di;
	double lat0=-43.56, lon0=172.12, dep0=0.0;
	int Nlat=20, Nlon=40, Ndep=10, NP=Nlat*Nlon*Ndep;
	double *latgrid, *longrid, *depgrid;
	int *allp;
	double dlat=1.2, dlon=1.2, ddep=4.0;
	time_t time1, time2;

	latgrid=dvector(1,NP);
	longrid=dvector(1,NP);
	depgrid=dvector(1,NP);
	allp=ivector(1,NP);
	crst.x=dvector(1,NP);
	crst.y=dvector(1,NP);

	double strikes[]={90, 90, 90, 90, 90, 90};
	double dips[]={90, 90, 45, 15, 15, 15};
	double rakes[]={0, 180, 45, 90, 270, 45};

	for (int i=1; i<=NP; i++) allp[i]=i;

	for (int i=1; i<=Nlat; i++){
		for (int j=1; j<=Nlon; j++){
			for (int k=1; k<=Ndep; k++){
				latgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lat0-0.5*dlat+(dlat/Nlat)*(i-0.5);
				longrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= lon0-0.5*dlon+(dlon/Nlon)*(j-0.5);
				depgrid[Nlat*Nlon*(k-1)+Nlat*(j-1)+i]= (ddep/Ndep)*(k-0.5);
			}
		}
	}

	init_crst(&crst);
	//broken: function does not exist anymore
	//read_farfalle_crust("/home/des/camcat/Code/CRS_2.01/input/inCan.dat", &crst);
	sprintf(fname0, "%s/okada/fm_slipmodels_summary.txt",testfolder);
	fout=fopen(fname0,"w");

	fprintf(fout, "strike\t dip \t rake \t magnitude \t resolution \n");
	eqfm.lat=lat0;
	eqfm.lon=lon0;
	eqfm.whichfm=1;
	eqfm.taper=ivector(1,4);
	for (int i=1; i<=3; i++) eqfm.taper[i]=1;
	eqfm.taper[4]=0;

	dcfs.NF=1;
	crst.nLat=Nlat;
	crst.nLon=Nlon;
	crst.nD=Ndep;
	crst.lat0=lat0;
	crst.lon0=lon0;
	crst.lat=latgrid;
	crst.lon=longrid;
	crst.depth=depgrid;
	crst.dlat=dlat/Nlat;
	crst.dlon=dlon/Nlon;
	crst.ddepth=ddep/Ndep;
	crst.x=dvector(1,NP);
	crst.y=dvector(1,NP);
	for (int i=1; i<=NP; i++) latlon2localcartesian(latgrid[i], longrid[i], lat0, lon0, crst.y+i, crst.x+i);
	latlon2localcartesian(eqfm.lat, eqfm.lon, crst.lat0, crst.lon0, &(eqfm.y), &(eqfm.x));

//	for (int i=1; i<=3; i++){
//		for (int j=1; j<=3; j++) crst.S[i][j]=0.0;
//	}

	time(&time1);
	for (int n=0; n<N; n++){
		eqfm.depth=dep0;
		eqfm.str1=strikes[n];
		eqfm.dip1=dips[n];
		eqfm.rake1=rakes[n];
		eqfm.whichfm=1;
		eqfm.mag=6.3;
		//res=0.5+ran1(&seed);

		fprintf(fout, "%.0lf\t%.0lf\t%.0lf\t%.1lf\t%.1lf\n", eqfm.str1, eqfm.dip1, eqfm.rake1, eqfm.mag, res);
		focmec2slipmodel(crst, &eqfm, res, 1, 1);
		find_gridpoints_d(crst.y, crst.x, depgrid, (int *) 0, 0, NP, eqfm.y, eqfm.x, eqfm.depth, eqfm.mag, 100000,  &(eqfm.nsel), &(eqfm.selpoints));
		dcfs.nsel=eqfm.nsel;
		dcfs.which_pts=eqfm.selpoints;
		printf("%d\n",dcfs.nsel);
		dcfs.S=d3tensor(1,dcfs.nsel, 1,3,1,3);
		dcfs.cmb=dvector(1,dcfs.nsel);
		sprintf(fname, "%s/okada/slipmodel%d.dat",testfolder, n);
		//print_slipmodel(fname, &eqfm, 1);

		// 1) using okadaDCFS:

		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;
		okadaDCFS(dcfs, &eqfm, 1, crst, &eqfm.str1, &eqfm.dip1, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, &(eqfm.rake1), 0);
		sprintf(fname, "%s/okada/okadaDCFSfixA%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

	//	okadaDCFS(dcfs, &eqfm, 1, crst, 0.0, 0.0, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, NULL, 1);
		sprintf(fname, "%s/okada/okadaDCFSfix_optrakeB%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, &(eqfm.rake1), 0);
		sprintf(fname, "%s/okada/okadaDCFSfixB%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;
		DCFScmbopt(&dcfs, 0, crst);
		sprintf(fname, "%s/okada/okadaDCFSopt%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		//2) using okadaCoeff:

		coeffs.NgridT=dcfs.nsel;
		coeffs.NP=eqfm.np_di*eqfm.np_st;
		okadaCoeff(&(coeffs.Coeffs_st), &(coeffs.Coeffs_dip), &eqfm, 1, crst, latgrid, longrid, depgrid);

		okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst,  &eqfm.str1, &eqfm.dip1, 0);
		sprintf(fname, "%s/okada/okadaCoeff_fix_optrakeA%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst, NULL, NULL, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, NULL, 1);
		sprintf(fname, "%s/okada/okadaCoeff_fix_optrakeB%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaCoeff2DCFS(coeffs.Coeffs_st, coeffs.Coeffs_dip, dcfs, &eqfm, crst, NULL, NULL, 1);
		resolve_DCFS(dcfs, crst, &eqfm.str1, &eqfm.dip1, &(eqfm.rake1), 0);
		sprintf(fname, "%s/okada/okadaCoeff_fix_A%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);
		for (int i=1; i<=dcfs.nsel; i++) dcfs.cmb[i]=0.0;

		okadaCoeff_resolve(coeffs, &resCoeff_st, &resCoeff_di, crst, &eqfm.str1, &eqfm.dip1, &eqfm.rake1);
		resolvedCoeff2DCFS(resCoeff_st, resCoeff_di, dcfs, &eqfm, crst);
		sprintf(fname, "%s/okada/okadaCoeff_fix_B%d.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);

		free_dvector(dcfs.cmb,1,dcfs.nsel);
		free_ivector(dcfs.which_pts,1,dcfs.nsel);
		free_d3tensor(dcfs.S,1,dcfs.nsel,1,3,1,3);
		free_f3tensor(coeffs.Coeffs_dip,1,coeffs.NP, 1, coeffs.NgridT, 1,6);
		free_f3tensor(coeffs.Coeffs_st,1,coeffs.NP, 1, coeffs.NgridT, 1,6);
		free_matrix(resCoeff_st, 1,coeffs.NP, 1, coeffs.NgridT);
		free_matrix(resCoeff_di, 1,coeffs.NP, 1, coeffs.NgridT);


	}
	fclose(fout);

	time(&time2);
	printf("Execution time: %.3lf seconds.\n", difftime(time2, time1));

	free_dvector(latgrid,1,NP);
	free_dvector(longrid,1,NP);
	free_dvector(depgrid,1,NP);
	free_dvector(crst.x,1,NP);
	free_dvector(crst.y,1,NP);
	free_ivector(allp,1,NP);

	return 0;

}

int test_latlon2localcartesian(){

	double n, e;

	latlon2localcartesian(0.0, 179.9, 0.0, -179.9,&n, &e);
	printf("-179.9 to 179.9: %lf\n",e);
	latlon2localcartesian(0.0, -179.9, 0.0, 179.9,&n, &e);
	printf("179.9 to -179.9: %lf\n",e);
	latlon2localcartesian(0.0, 179.9, 0.0,  180.1,&n, &e);
	printf("181.1 to 179.9: %lf\n",e);
	latlon2localcartesian(0.0, -179.9, 0.0, -180.1,&n, &e);
	printf("-181.1 to -179.9: %lf\n",e);

	return (0);

}

int test_distance(){

	 /* Results can be visualized in matlab:
	 *  for n=1:4
		figure(n)
		Plot_sliponfault(n,strcat('~/Code/dist2fault/tests/',num2str(n),'_fault_'))
		l=load(strcat('~/Code/CRS_2.01/test/',num2str(n),'_dist.dat'));
		scatter3(l(:,2),l(:,1),-l(:,3),30,abs(l(:,4)),'filled')
		caxis([-100 100])
		end
	 */

	FILE *fout1;
	char fname1[120];
	long seed=-36294638;
	double strikes[]={0.0, 45.0, 0.0, 0.0};
	double dips[]={90.0, 70.0, 15.0, 0.1};
	double Lat0=45, Lon0=0.0, D0=0;
	double pos_d[]={0,50};
	double pos_s[]={-60,60};
	double *d;
	double dlat=2.0, dlon=2.0, ddep=100;
	int NP=1000;
	double lats[NP+1], lons[NP+1], deps[NP+1];

	for (int i=1; i<=NP; i++){
		lats[i]= Lat0-0.5*dlat+ dlat*ran1(&seed);
		seed=-seed;
		lons[i]= Lon0-0.5*dlon+ dlon*ran1(&seed);
		seed=-seed;
		deps[i]= ddep*ran1(&seed);
		seed=-seed;
	}


	for (int i=0; i<4; i++){
		d=dist2fault0(lats, lons, deps, NP, strikes[i], dips[i],	Lat0, Lon0, D0, pos_s, pos_d);

		sprintf(fname1, "%s/%d_fault_allpatches.dat", testfolder, i+1);
		fout1=fopen(fname1, "w");
		fprintf(fout1,  "%d   %.4lf   %.4lf   %.3lf   %.2lf   %.2lf   %.3lf   %.3lf   %d   %d   %.5lf\n", 1, Lat0, Lon0, D0, pos_s[1]-pos_s[0], pos_d[1]-pos_d[0], strikes[i], dips[i], 1, 1, 0.0);
		fclose(fout1);
		sprintf(fname1, "%s/%d_fault_patch1.dat", testfolder, i+1);
		fout1=fopen(fname1, "w");
		fprintf(fout1, "%12.5lf\t%12.5lf\t%12.5lf\t%12.5lf\t%12.5lf\n", 0.5*(pos_s[0]+pos_s[1]), 0.5*(pos_d[0]+pos_d[1]), 1.0, 1.0, 0.0);
		fclose(fout1);

		sprintf(fname1, "%s/%d_dist.dat", testfolder, i+1);
		fout1=fopen(fname1, "w");
		for (int p=1; p<=NP; p++) fprintf(fout1, "%lf\t%lf\t%lf\t%lf\n", lats[p], lons[p], deps[p], d[p]);
		fclose(fout1);

		free_dvector(d,1,NP);
	}

	printf("done.\n");
	return 0;
}

int test_matrix(){

//	double 	s[]={-10.0,-0.5,-0.1},\
//			st[]={115.0, 0.0, 25.0},\
//			di[]={0.0,90.0,0.0};

	double 	s[]={-10.0,-0.1,-0.5},\
			st[]={115.0, 25.0, 205.0},\
			di[]={0.0,10.0,80.0};

//	double  s[]={-5.0, 5.0, 0.0},\
//			st[]= {6.646198, -83.345443, 60.000013 },	\
//			di[]={-0.596911, -0.802279, 89.000000};

	float eig[4];
	float **v;
	double **S, *sigma;
	float **Sf;
	int j;
	double cmb, st1, st2, di1, di2, ra1, ra2;

	sigma=dvector(0,2);
	v=matrix(1,3,1,3);
	Sf=matrix(1,3,1,3);
	S=prestress_eigen(s, st, di);
	printf("\nS: \n");
	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++) {
			printf("%lf\t", S[i][j]);
			Sf[i][j]=(float)S[i][j];
		}
		printf("\n");
	}

	jacobi(Sf, 3, eig, v, &j);
	printf("\nv: \n");
	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++) {
			printf("%lf\t", v[i][j]);
		}
		printf("\n");
	}

	printf("\n eigenvalues:\n");
	for (int i=1; i<=3; i++) printf("%f\n", eig[i]);


	//transform eigenvectors into strike, dip:
	for (int i=1; i<=3; i++){
		di[i-1]=RAD2DEG*asin(v[3][i]);
		st[i-1]=RAD2DEG*atan2(v[2][i],v[1][i]);
		sigma[i-1]=(double) eig[i];
	}

	//st[2]=0.0;

	printf("\nstrikes: %lf, %lf, %lf \n", st[0], st[1], st[2]);
	printf("\ndips:    %lf, %lf, %lf \n", di[0], di[1], di[2]);

	//recalculate S from eigenvectors:
	S=prestress_eigen(sigma, st, di);
	printf("\nS: \n");
	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++) {
			printf("%lf\t", S[i][j]);
			Sf[i][j]=(float)S[i][j];
		}
		printf("\n");
	}


	cmbopt(S[1][1], S[2][2], S[3][3], S[1][2], S[2][3], S[1][3], 0.0, 0.4, 0.0, 90.0, 180.0, &cmb, &st1, &di1, &ra1, &st2, &di2, &ra2);
	printf("\n oops:\n");
	printf("st1=%lf\t di1=%lf\t ra1=%lf\n", st1, di1, ra1);
	printf("st2=%lf\t di2=%lf\t ra2=%lf\n", st2, di2, ra2);

	return 0;
}

int test_matrix2(){

	double 	s[]={5.0,-5.0,0.0};
	double st[3], di[3], sigma[3];
	double st_oop=330.0;
	double di_oop=89.0;
	double ra_oop=180;
	double p=0.0, fr=0.3;

	float eig[4];
	float **v;
	double **S;
	float **Sf;
	int j;
	double cmb, st1, st2, di1, di2, ra1, ra2;

	v=matrix(1,3,1,3);
	Sf=matrix(1,3,1,3);

	//calculate matrix using orientation of oops:
	prestress(s[0],s[1],s[2],st_oop, di_oop, ra_oop, p, fr, &S);
	printf("\nS: \n");
	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++) {
			printf("%lf\t", S[i][j]);
			Sf[i][j]=(float)S[i][j];
		}
		printf("\n");
	}

	//find eigenvalues/eigenvectors of matrix:
	jacobi(Sf, 3, eig, v, &j);
	printf("\nv: \n");
	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++) {
			printf("%lf\t", v[i][j]);
		}
		printf("\n");
	}

	printf("\n eigenvalues:\n");
	for (int i=1; i<=3; i++) printf("%f\n", eig[i]);

	//transform eigenvectors into strike, dip:
	for (int i=1; i<=3; i++){
		di[i-1]=RAD2DEG*asin(v[3][i]);
		st[i-1]=RAD2DEG*atan2(v[2][i],v[1][i]);
		sigma[i-1]=(double) eig[i];
	}

	printf("\nstrikes: %lf, %lf, %lf \n", st[0], st[1], st[2]);
	printf("\ndips:    %lf, %lf, %lf \n", di[0], di[1], di[2]);

	//recalculate S from eigenvectors:
	S=prestress_eigen(sigma, st, di);
	printf("\nS2: \n");
	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++) {
			printf("%lf\t", S[i][j]);
			Sf[i][j]=(float)S[i][j];
		}
		printf("\n");
	}

	//calculate oops from new S:
	cmbopt(S[1][1], S[2][2], S[3][3], S[1][2], S[2][3], S[1][3], 0.0, fr, 0.0, 90.0, 180.0, &cmb, &st1, &di1, &ra1, &st2, &di2, &ra2);
	printf("\n oops:\n");
	printf("st1=%lf\t di1=%lf\t ra1=%lf\n", st1, di1, ra1);
	printf("st2=%lf\t di2=%lf\t ra2=%lf\n", st2, di2, ra2);

	return 0;
}

int test_readZMAP(){

	char xmlfile[]="input/nz-forecast-template-M4.xml";
	char file[]="input/nz_2009_2011_fixed.dat";
	char fore_file[]="input/darf_temp.txt";
	char crust_file[]="input/inCan.dat";
	//char file[]="input/quake_recent.dat";
	char output[120];
	struct catalog cat;
	struct eqkfm *eqfm;
	struct crust crst;
	struct tm reftime, tnow;
	FILE *fout;

	double t0s=-10000, t1s=10000, t0c=-10000, t1c=10000;
	double border=5;
	double dDCFS=100;
	int junk;
	double djunk;
	double res=3.0, res_z=1.0;

	sscanf("2010-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	reftime.tm_year-=1900;
	reftime.tm_mon-=1;
	reftime.tm_isdst=0;

//	read_xmltemplate(xmlfile, &tnow, &t0c, &t1c, &junk, &junk, &crst, &djunk, &djunk);
	//crst.x=crst.y=crst.dAgrid=crst.depth;

	//broken (no crust_file anymore)
	//read_crust(crust_file, fore_file, NULL, &crst, res, res_z);
	readZMAP (&cat, &eqfm, NULL, file, crst, reftime, t0s, t1s, t0c, t1c, 10.0, 0.0, border, border, dDCFS, 0);

//	readZMAP (&cat, &eqfm, NULL, file, crst, tnow, t0s, t1s, t0c, t1c, 10.0, 0.0, border, border, dDCFS, 1);
//	sprintf(output,"%s/Darfield_cat.dat", testfolder);
//	fout=fopen(output,"w");
//	for (int i=1; i<=cat.Z; i++) fprintf(fout,"%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", i, cat.t[i], cat.lat0[i], cat.lon0[i], cat.depths0[i], cat.mag[i]);
//	fclose(fout);

	return (0);
}

int test_assign_GRnorm(){

	int N=6;
	double mags[]={1.5,2.5,3.5,4.5,5.5,6.5};
	double b=0.8;
	double *res;

	res=assign_GRnorm(mags, N, b, 1);

	for (int i=0; i<N; i++) printf("%.1lf\t%.6lf\n",mags[i], res[i]);

	return 0;

}

int test_bin_equnumber(){

	int N=20, Nbin=10;
	double v[]={0.5, 0.7183, 0.9446, 0.0878, 0.2795, 0.5962, 0.8284, 0.7822, 0.5572, 0.0363, 0.6694, 0.8490, 0.0655, 0.3608, 0.2581, 0.4326, 0.3061, 0.9666, 0.1299, 0.2174};
	double *bc, *bh;

	bin_equnumber(v, N, Nbin, &bc, &bh);

	return 0;
}

int test_Mc_maxcurv(){

	char inputfile[]="input/nz_merged.dat";
	int c=8;
	long r;
	double **data;
	double *mags;
	double Mc, b;

	r=countline(inputfile)+1;

	data=dmatrix(1,c,1,r);
	mags=data[6];
	read_matrix(inputfile, c, 0, data, &r);

	Mc=Mc_maxcurv(mags, (int)r);
	b= calculatebvalue(mags, (int) r, Mc);

	printf("Mc=%.3lf\n b=%.3lf\n", Mc,b);
	return 0;
}

//int test_xml(){
//
////	char *value;
////	fetch_xml_value("input/nz_temp.forecast.xml", "depthLayer", &value);
////	printf("%s\n", value);
//
//	char **values;
//	int N, r;
//	struct crust crst;
//
//	//r=fetch_xml_property("input/nz_temp.forecast.xml", "depthLayer", "max", &values, &N);
//	//r=fetch_xml_property2("input/nz_temp.forecast.xml", "cell", "lon", &values, &N);
//	//r=fetch_xml_property3("input/nz_temp.forecast.xml", "bin", "m", &values, &N);
////	for (int i=0; i<N; i++) printf("value=%s\n", values[i]);
////	printf("return value: %d (N=%d)\n", r, N);
//
//	double tn, t0, t1, m0, m1;
//	double *rates, *GR;
//	int Minf, nobins;
//	char *inputfile="input/nz_temp.forecast.xml";
//	char *outputfile="output/nz_temp.forecast.xml";
//
//	read_xmltemplate(inputfile, NULL, NULL, NULL, &nobins, &Minf, &crst, &m0, &m1);
//
//	rates=dvector(1,crst.N_allP);
//	GR=dvector(1,nobins);
//
//	for (int i=1; i<=crst.N_allP; i++) rates[i]=10.0*i;
//	for (int i=1; i<=nobins; i++) GR[i]=(double)i;
//
//	xml_writeforecast(inputfile, outputfile, rates, GR, "%.1lf");
//
//	return;
//}

void test_convertgeometry(){
//test ok.
//update, 29.10: also test conversion from low to high resolution. test ok.

	struct crust cr, cr1;
	struct pscmp *d, *d1;
	double *o, *n=0, *n1=0;
	int err, ind;
	int P[3], Pn[3], nsub[3];
	char fname[120];
	extern int farfalle;

//	read_crust("input/crust.dat", NULL, &cr, 3.0,5.0);

	//broken (no crust_file anymore)
	//read_crust("input/inCan.dat", "input/darf_temp.txt", NULL, &cr, 3.0,5.0);
	d=pscmp_arrayinit(cr, 0,0);	//used for output.
	d[0].which_pts=ivector(1,cr.N_allP);

	o=dvector(1,cr.N_allP);
	n1=dvector(1,cr.N_allP);

	nsub[0]=cr.nLat/cr.nLat_out;
	nsub[1]=cr.nLon/cr.nLon_out;
	nsub[2]=cr.nD/cr.nD_out;

	//assign value given by coordinate that will be in new geometry (easier to check):
	for (int i=1; i<=cr.N_allP; i++) {
		d[0].which_pts[i]=i;
		P[0]=((i-1)%(cr.nLat*cr.nLon))%cr.nLat+1;
		P[1]=(((i-1)%(cr.nLat*cr.nLon))+1-P[0])/cr.nLat +1;
		P[2]=(i-P[0]-(P[1]-1)*cr.nLat)/(cr.nLat*cr.nLon) +1;
		for (int n=0; n<3; n++) Pn[n]=(P[n]-1)/nsub[n]+1;
		o[i]=(double) (Pn[0]+cr.nLat_out*(Pn[1]-1)+(cr.nLat_out*cr.nLon_out)*(Pn[2]-1));
//		printf("[%d, %d, %d] \t %.0lf\n",Pn[0], Pn[1], Pn[2], o[i]);
	}

	err=convert_geometry(cr, o, &n, 0, 0);

	sprintf(fname, "%s/old_grid2910.dat", testfolder);
	d[0].nsel=cr.N_allP;
	print_grid(fname,d[0],cr, o);

	//copy new geometry into a pscmp structure: warning! may be wrong (change quickly after removing lat, lon depth from pscmp structure).
	cr1.nD=cr.nD_out;
	cr1.nLat=cr.nLat_out;
	cr1.nLon=cr.nLon_out;
	cr1.N_allP=cr.nLat_out*cr.nLon_out*cr.nD_out;
	cr1.lat=dvector(1,cr1.N_allP);
	cr1.lon=dvector(1,cr1.N_allP);
	cr1.depth=dvector(1,cr1.N_allP);
	cr1.dlat=cr.dlat_out;
	cr1.dlon=cr.dlon_out;
	cr1.ddepth=cr.ddepth_out;

	for (int dd=1; dd<=cr1.nD; dd++){
		for (int lo=1; lo<=cr1.nLon; lo++){
			for (int la=1; la<=cr1.nLat; la++){
				ind=(dd-1)*cr1.nLon+(lo-1)*cr1.nLat+la;
				cr1.lat[ind]=cr.latmin+(la-0.5)*cr1.dlat;
				cr1.lon[ind]=cr.lonmin+(lo-0.5)*cr1.dlon;
				cr1.depth[ind]=cr.depmin+(dd-0.5)*cr1.ddepth;
			}
		}
	}

	d1=pscmp_arrayinit(cr1, 0,0);
	d1[0].nsel=cr1.N_allP;
	d1[0].which_pts=ivector(1,cr1.N_allP);
	for (int i=1; i<=cr1.N_allP; i++) d1[0].which_pts[i]=i;

	sprintf(fname, "%s/new_grid2910.dat", testfolder);
	print_grid(fname, d1[0],cr1, n);

	err=convert_geometry(cr, n, &n1, 1, 1);
	sprintf(fname, "%s/old_grid_rec2910.dat", testfolder);
	print_grid(fname, d[0],cr, n1);

	printf("Test completed!\n");

	return;
}

//void test_readslipmodelold(){
//
//	struct eqkfm *eq;
//	int NF;
//
//	eq=eqkfm_array(0,10);
//
//	read_justslipmodel("input/slip_dar_beav_pscmp.dat", eq, &NF, 3000);
//}

int testspeed_coeff(){
	struct eqkfm eqkfm0;
	struct Coeff_LinkList AllCoeff;
	struct crust crst;
	char slipmodel[120];
	int nlat, nlon, ndep, ind, NF;
	double lat0, lon0, dep0, Dlat=2.5, Dlon=2.5, Ddep=-40, dlat=0.1, dlon=0.1, ddep=1.0;
	double *lats, *lons, *deps;
	time_t 	extime0, extime1;
	FILE *f1, *f2;

	crst.lambda=31226, crst.mu=26624;
	sprintf(slipmodel,"/home/des/camcat/Code/CRS_1.0/INPUT/Parkfield_new.inp");
	/*read_justslipmodel(slipmodel, &eqkfm0, &NF, crst.mu);*/		//obsolete function -> test doesn't work anymore.

	lat0=eqkfm0.lat-0.5*Dlat+0.5*dlat;
	lon0=eqkfm0.lon-0.5*Dlon+0.5*dlon;
	dep0=0.5*ddep;

	nlat=(int) (Dlat/dlat);
	nlon=(int) (Dlon/dlon);
	ndep=(int) -(Ddep/ddep);

	eqkfm0.nsel=nlat*nlon*ndep;
	eqkfm0.selpoints=ivector(1,eqkfm0.nsel);
	for (int i=1; i<=eqkfm0.nsel; i++) eqkfm0.selpoints[i]=i;


	lats=dvector(1,nlat*nlon*ndep);
	lons=dvector(1,nlat*nlon*ndep);
	deps=dvector(1,nlat*nlon*ndep);

	for (int d=1; d<=ndep; d++)	{
		for (int lo=1; lo<=nlon; lo++){
			for (int la=1; la<=nlat; la++){
				ind=(d-1)*nlat*nlon+(lo-1)*nlat+la;
				lats[ind]=lat0+(la-1.0)*dlat;
				lons[ind]=lon0+(lo-1.0)*dlon;
				deps[ind]=dep0+(d-1.0)*ddep;
			}
		}
	}

	time(&extime0);
   	okadaCoeff(&(AllCoeff.Coeffs_st), &(AllCoeff.Coeffs_dip), &eqkfm0, 1, crst, lats, lons, deps);
	time(&extime1);
	printf("Time to calculate coefficients: %f sec\n", difftime(extime1,extime0));

	time(&extime0);
	f1=fopen("tests/Coeff_st.dat","w");
	f2=fopen("tests/Coeff_di.dat","w");

	for (int p=1; p<=eqkfm0.np_di*eqkfm0.np_st; p++){
		for (int n=1; n<=eqkfm0.nsel; n++){
			for (int i=1; i<=6; i++){
				fprintf(f1,"%.5f\t", AllCoeff.Coeffs_st[p][n][i]);
				fprintf(f2,"%.5f\t", AllCoeff.Coeffs_dip[p][n][i]);
			}
		}
	fprintf(f1,"\n");
	fprintf(f2,"\n");
	}
	fclose(f1);
	fclose(f2);
	time(&extime1);
	printf("Time to write to file: %f sec\n", difftime(extime1,extime0));

	time(&extime0);
	f1=fopen("tests/Coeff_st.dat","r");
	f2=fopen("tests/Coeff_di.dat","r");

	for (int p=1; p<=eqkfm0.np_di*eqkfm0.np_st; p++){
		for (int n=1; n<=eqkfm0.nsel; n++){
			for (int i=1; i<=6; i++){
				fscanf(f1,"%f\t", &(AllCoeff.Coeffs_st[p][n][i]));
				fscanf(f2,"%f\t", &(AllCoeff.Coeffs_dip[p][n][i]));
			}
		}
	}


	fclose(f1);
	fclose(f2);
	time(&extime1);
	printf("Time to read from file: %f sec\n", difftime(extime1,extime0));

   	return(0);
   	}

int test_hash(){

	long int res2, res=1;
	char string[120]="blablablablablabla";
	char string2[120]="blablablablablablablu";

	sprintf(string2,"%s%s",string,string);
	res=hashlittle( string, strlen(string), 1);
	res2=hashlittle( string2, strlen(string), 1);

	printf("%ld\n",res);
	printf("%ld\n",res2);

	printf("\n string length= %d\n", (int) strlen(string));
	printf("\n string2 length= %d\n", (int) strlen(string2));

	return (res==res2);
}

void test_taper_multislip(){
/*
 */
	int N=2;
	long seed=-19329935;
	double res=1.5, mag=7.1, sign, sign2;
	double north, east;
	struct eqkfm eqfm[3];
	struct crust crst;
	char fname[120];
	double d_min=3.0;
	double noise1, noise2;

	crst.lambda=31226, crst.mu=26624;//calculated for Vp=5.7,Vs=3.2, rho=2600 (from Wang psgrn input file for Parkfield). MPa.

	for (int n=1; n<N; n++){
		sprintf(fname, "%s/pseudo_Darfield_%d.dat",testfolder, n+N-1);
		eqfm[0].lat=-43.56;
		eqfm[0].lon=172.12;
		eqfm[0].depth=0.0;
		eqfm[0].whichfm=1;
		eqfm[0].str1=95.0-5.0+10.0*ran1(&seed);
		eqfm[0].dip1=90.0-2.0+4.0*ran1(&seed);
		eqfm[0].rake1=180.0;
		eqfm[0].mag=mag;
		if (n==1) eqfm[0].taper=ivector(1,4);
		for (int i=1; i<=4; i++) eqfm[0].taper[i]=0;

		focmec2slipmodel(crst, eqfm, res, 1, 1);
		sign2=ran1(&seed)-0.5;
		sign2*=1.0/fabs(sign2);	//-1 or 1;

		for (int f=1; f<3; f++){
			eqfm[f].str1=95.0-5.0+10.0*ran1(&seed);
			eqfm[f].dip1=90.0-2.0+4.0*ran1(&seed);
			eqfm[f].rake1=180.0;
			eqfm[f].whichfm=1;
			eqfm[f].mag=mag;
			if (n==1) eqfm[f].taper=ivector(1,4);
			for (int i=1; i<=4; i++) eqfm[f].taper[i]=0;
			sign=sign2;
			sign2=ran1(&seed)-0.5;
			sign2*=1.0/fabs(sign2);	//-1 or 1;

			noise1= 0.0*ran1(&seed);
			noise2= 0.0*ran1(&seed);

			focmec2slipmodel(crst, eqfm+f, res, 1, 1);
			north= 0.5*(sign*eqfm[f-1].L*cos(DEG2RAD*eqfm[f-1].str1)+sign2*eqfm[f].L*cos(DEG2RAD*eqfm[f].str1))+noise1;
			east= 0.5*(sign*eqfm[f-1].L*sin(DEG2RAD*eqfm[f-1].str1)+sign2*eqfm[f].L*sin(DEG2RAD*eqfm[f].str1))+noise2;
			localcartesian2latlon(north, east, eqfm[f-1].lat, eqfm[f-1].lon,  &(eqfm[f].lat), &(eqfm[f].lon));
			//focmec2slipmodel(crst, eqfm+f, res, 1, 1);

		}

		//which_taper(eqfm,  3, 1, 1, d_min);
		//for (int f=0; f<3; f++) suomod1_taper(eqfm[f], eqfm+f);

		print_slipmodel(fname, eqfm, 3);
	}

}

void test_focmec2slipmodel(){
/* to visualize results, should then run the following:
 * for k in $(seq 0 9); do  /home/des/camcat/Code/Scripts/slipmodel_sep.sh fm_slipmodel_$k.dat ; done
 * and matlabl script test_focmec2slipmodel.
 */

	int N=10;
	long seed=-19329935;
	double res;
	struct eqkfm eqfm;
	struct crust crst;
	FILE *fout;
	char fname[120], fname0[120];

	crst.lambda=31226, crst.mu=26624;//calculated for Vp=5.7,Vs=3.2, rho=2600 (from Wang psgrn input file for Parkfield). MPa.

//	sprintf(fname0, "%s/fm_slipmodels_summary.txt",testfolder);
//	fout=fopen(fname0,"w");

//	fprintf(fout, "strike\t dip \t rake \t magnitude \t resolution \n");
	eqfm.lat=-43.56;
	eqfm.lon=172.12;
	eqfm.whichfm=1;
	eqfm.taper=ivector(1,4);
	for (int i=1; i<=3; i++) eqfm.taper[i]=1;
	eqfm.taper[4]=0;

	for (int n=0; n<N; n++){
		eqfm.str1=360*ran1(&seed);
		eqfm.dip1=90*ran1(&seed);
		eqfm.rake1=360*ran1(&seed);
		eqfm.whichfm=1;
		eqfm.mag=6.0+2.0*ran1(&seed);
		//res=0.5+ran1(&seed);
		res=3.0;


//		fprintf(fout, "%.0lf\t%.0lf\t%.0lf\t%.1lf\t%.1lf\n", eqfm.str1, eqfm.dip1, eqfm.rake1, eqfm.mag, res);
		focmec2slipmodel(crst, &eqfm, res, 1, 1);

		sprintf(fname, "%s/suomod_minus_%d.dat",testfolder, n);
		print_slipmodel(fname, &eqfm, 1);
	}
//	fclose(fout);

}

void tests_eqkfm_addslipmodels(){

	struct slipmodels_list all_slipmodels;
	struct eqkfm *eq_in1, *eq_out;
	struct crust crst;
	int N1=10, Nout;
	int *indices;
	double tmain2[]={2.0, 3.0, 7.5, 8.0, 9.0};
	double mmain2[]={4.4, 4.6, 4.4, 4.6, 4.8};
	//double mmain2={2.5, 3.0, 6.0, 7.3, 8.0};
	int Nfaults2[]={3,3,3,3,3};
	int pts[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
	int err;
	int *NFout;

	crst.lambda=31226, crst.mu=26624;//calculated for Vp=5.7,Vs=3.2, rho=2600 (from Wang psgrn input file for Parkfield). MPa.
	crst.N_allP=20;
	crst.list_allP=pts-1;


	all_slipmodels.is_afterslip=0;
	all_slipmodels.tmain=tmain2;
	all_slipmodels.mmain=mmain2;
	all_slipmodels.Nfaults=Nfaults2;
	all_slipmodels.NSM=5;
	all_slipmodels.no_slipmodels=ivector(0,all_slipmodels.NSM-1);
	all_slipmodels.disc=dvector(0,all_slipmodels.NSM-1);
	all_slipmodels.slipmodels=	malloc((all_slipmodels.NSM-1) * sizeof(char*));
	for (int nn=0; nn<all_slipmodels.NSM; nn++) {
		all_slipmodels.slipmodels[nn] = malloc(120 * sizeof(char));
	}

	eq_in1=eqkfm_array(0,N1-1);
	for (int i=0; i<N1; i++){
		eq_in1[i].taper=ivector(1,4);
		for (int ii=1; ii<=4; ii++) eq_in1[i].taper[ii]=1;
		eq_in1[i].str1=360*(double) (i%N1);
		eq_in1[i].dip1=90*(double) (i+5%N1);
		eq_in1[i].rake1=180.0;
		eq_in1[i].whichfm=1;
		eq_in1[i].t=(double) i;
		eq_in1[i].mag=4.0+0.2*(double)(i%5);		//[4.0, 4.2, 4.4, 4.6, 4.8, 4.0, ...];
		eq_in1[i].is_mainshock= (eq_in1[i].mag>=4.3);
	}

	for (int i=0; i<all_slipmodels.NSM; i++){

		for (int c=0; c<all_slipmodels.Nfaults[i]; c++){
			all_slipmodels.no_slipmodels[i]=1;
			all_slipmodels.disc[i]=1.0;
			sprintf(all_slipmodels.slipmodels[i],"input/pseudo_Darfield/pseudo_Darfield_1.dat");
		}
	}

	err= eqkfm_addslipmodels(eq_in1, all_slipmodels, &eq_out, &indices, N1, &Nout, &NFout, 0.001, 0.3, 1.0, crst, 1, 1);

	return;
}

void test_suomod1_hf(){
	/* to visualize (in Matlab):
	 *
	 for k=1:9; figure(k); Plot_sliponfault(k,strcat('~/Code/CRS_2.0/tests/suomod1_',num2str(k),'_')); set(gca,'CameraPosition',[137.4177  -43.6711   -9.0672]);  end
	 *
	 */

	struct eqkfm eqfm, eqfm2;
	struct crust crst;
	long seed=-19329935;
	char fname0[120];
	double H=1.0;	//Hurst exponent.
	int N=10;
	int noise_only=0;
	time_t t0,t1;

	crst.lambda=31226, crst.mu=26624;//calculated for Vp=5.7,Vs=3.2, rho=2600 (from Wang psgrn input file for Parkfield). MPa.

	eqfm.lat=-43.56;
	eqfm.lon=172.12;
	eqfm.whichfm=1;
	eqfm.taper=ivector(1,4);
	for (int i=1; i<=4; i++) eqfm.taper[i]=1;
	eqfm.whichfm=1;
	eqfm.mag=7.0;
	eqfm.str1=0.0;
	eqfm.dip1=90.0;
	eqfm.rake1=180.0;

	focmec2slipmodel(crst, &eqfm, 1.0,1,0);

//	time(&t0);
//	for (int n=1; n<=N; n++){
//		sprintf(fname0, "%s/suomod1old_%d.dat",testfolder,n);
//		suomod1_addhf_old(eqfm, &eqfm2, 0.0, &seed, 1,1);
//		print_slipmodel(fname0, &eqfm2, 1);
//	}
//	time(&t1);
//	printf("Old function takes %.2f seconds\n", difftime(t1,t0));

	time(&t0);
	for (int n=1; n<=N; n++){
		sprintf(fname0, "%s/suomod%.1f_%d.dat",testfolder,H,n);
		suomod1_hf(eqfm, &eqfm2, H, &seed, noise_only);
		print_slipmodel(fname0, &eqfm2, 1);
	}
	time(&t1);

	//printf("New function takes %.2f seconds\n", difftime(t1,t0));

}

void test_reduction(){
	int   i, n, chunk;
	float a[100], b[100], result;

	/* Some initializations */
	n = 100;
	chunk = 10;
	result = 0.0;
	for (i=0; i < n; i++)
	{
		a[i] = i * 1.0;
		b[i] = i * 2.0;
	}

	#pragma omp parallel for      \
	default(shared) private(i)  \
	schedule(static,chunk)      \
	reduction(+:result)

	for (i=0; i < n; i++) result = result + (a[i] * b[i]);

	printf("Final result= %f\n",result);
	return;
}

void make_random_catalog(){

	long seed=-18483024;
	char *fname="input/dummy_catalog.dat";
	FILE *fout;
	int N=100;
	double 	t0=0.0,\
			t1=20.0,\
			lat1=-44.76,\
			lat2=-42.34,\
			lon1=171.04,\
			lon2=173.46,\
			d1=0,\
			d2=44,\
			m1=3.0,\
			m2=7.0,\
			sdh=1.0,\
			sdd=3.0;
	double lat, lon, d, m, t, dt=(t1-t0)/(1.0*N);

	fout=fopen(fname,"w");
	for (int i=1; i<=N; i++){
		lat=lat1+ran1(&seed)*(lat2-lat1);
		lon=lon1+ran1(&seed)*(lon2-lon1);
		d=d1+ran1(&seed)*(d2-d1);
		m=m1+ran1(&seed)*(m2-m1);
		t=t0+(i-1)*dt;
		fprintf(fout, "%d\t%.6lf\t%.3lf\t%.3lf\t%.3lf\t%.2lf\t%.1lf\t%.1lf\n", i, t, lat, lon, d, m, sdh, sdd);
	}
	fclose(fout);
	return;
}
