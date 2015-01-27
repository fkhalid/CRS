/*
 * background_rates.c
 *
 *  Created on: Nov 26, 2013
 *      Author: camcat
 */

#include "d_background_rates.h"

void d_background_rates(){
/*  Produces background rates varying the time period, to test how long a catalog should be used.
 */
	//char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/MaxWerner.Canterbury_M2_1980_2012.txt";
	char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/1900_2011.9.3.dat";
	char fname[120];
	FILE *fout;
	char crust_file[]="input/inCan.dat";
	char fore_file[]="input/darf_temp.txt";
	struct crust crst;
	struct tm reftime;
	int Nt=10;
	double t0, t00=100*365;
	double Mcut=20.0, Mmain=6.5;
	double dR=50, dZ=50;
	int ord=2, err;
	double smoothing=3.0;	//min. distance used for smoothing.
	double res=3.0, res_z=1.0;

	sscanf("2010-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	reftime.tm_year-=1900;
	reftime.tm_mon-=1;
	reftime.tm_isdst=0;

	read_crust(crust_file, fore_file, &crst, res, res_z);

	fout=fopen("test/bg_rate_100yrsS3.dat","w");
	for (int t=1; t<=Nt; t++){
		t0=-t00+(t-1)*(t00/Nt);
		err=background_rate2(cat_file, &crst, reftime, Mcut, Mmain, t0, 0.0, dR, dZ, smoothing, ord);
		sprintf(fname,"test/bg_rate_100yrsS3%.0lf.dat",t0);
		//print_rate(fname, crst, NULL);
		fprintf(fout,"%.3lf\t%.3lf\n",t0,crst.r0);
	}

	return;
}

void background_rates(){
/*  Produces background rates varying completeness magnitide.
 *
 */

	//char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/MaxWerner.Canterbury_M2_1980_2012.txt";
	//char cat_file[]="/home/des/camcat/Data/Catalogs/Others/jma_cat_2010_2013_update20130329.dat";
	char cat_file[]="/home/des/camcat/Data/Catalogs/Others/jma_cat_2010_update20130329.dat";
	char fname[120];
	char fore_file[]="input/other/tohoku_bgrate_nonuniform_3to39.dat";
	struct crust crst;
	struct tm reftime;
	//double t0=-434.614525;	//from 01/01/2010 (start of catalog).
	//double t0=-79.0;	//from 01/01/2011 (start of catalog).
	double Mmain=6.5;
	double dR=50, dZ=50;
	int ord=2, err=0;
	double smoothing=5.0;	//min. distance used for smoothing.
	double res=3.0, res_z=1.0;
	double Mcuts[6]={-1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
	int Nmcuts=6;

	//sscanf("2010-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	sscanf("2011-03-11T14:46:18Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	reftime.tm_year-=1900;
	reftime.tm_mon-=1;
	reftime.tm_isdst=0;

	err=read_crust(fore_file, NULL, &crst, res, res_z, 0);
	if (err) return;

	for (int i=0; i<Nmcuts; i++){
		sprintf(fname,"test/bgrate_toho2010_%.2f.dat",Mcuts[i]);
		err=background_rate2(cat_file, &crst, reftime, Mcuts[i], Mmain, NULL, NULL, 0, dR, dZ, smoothing, ord);
		//print_rate(fname, crst, NULL);
	}

	return;
}

