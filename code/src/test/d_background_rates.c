
/*   Copyright (C) 2015 by Camilla Cattania and Fahad Khalid.
 *
 *   This file is part of CRS.
 *
 *   CRS is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CRS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CRS.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 * background_rates.c
 *
 *  Created on: Nov 26, 2013
 *      Author: camcat
 */

#include "d_background_rates.h"

#include <stdio.h>
#include <time.h>

#include "../inp_out/write_csep_forecast.h"
//#include "../defines.h"
//#include "../inp_out/read_crust.h"
//#include "../seis/background_rate.h"

//void d_background_rates(){
///*  Produces background rates varying the time period, to test how long a catalog should be used.
// */
//	//char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/MaxWerner.Canterbury_M2_1980_2012.txt";
//	char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/1900_2011.9.3.dat";
//	char fname[120];
//	FILE *fout;
//	char crust_file[]="input/inCan.dat";
//	char fore_file[]="input/darf_temp.txt";
//	struct crust crst;
//	struct tm reftime;
//	int Nt=10;
//	double t0, t00=100*365;
//	double Mcut=20.0, Mmain=6.5;
//	double dR=50, dZ=50;
//	int ord=2, err;
//	double smoothing=3.0;	//min. distance used for smoothing.
//	double res=3.0, res_z=1.0;
//
//	sscanf("2010-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
//	reftime.tm_year-=1900;
//	reftime.tm_mon-=1;
//	reftime.tm_isdst=0;
//
//	read_crust(crust_file, fore_file, &crst, res, res_z);
//
//	fout=fopen("test/bg_rate_100yrsS3.dat","w");
//	for (int t=1; t<=Nt; t++){
//		t0=-t00+(t-1)*(t00/Nt);
//		err=background_rate(cat_file, &crst, reftime, Mcut, Mmain, t0, 0.0, dR, dZ, smoothing, ord);
//		sprintf(fname,"test/bg_rate_100yrsS3%.0lf.dat",t0);
//		//print_rate(fname, crst, NULL);
//		//fprintf(fout,"%.3lf\t%.3lf\n",t0,crst.r0);
//	}
//
//	return;
//}

void background_rates(){
/*  Produces background rates
 *
 */

	//char cat_file[]="/home/des/camcat/Data/Catalogs/ZMAP/new_zeland/MaxWerner.Canterbury_M2_1980_2012.txt";
	char cat_file[]="/home/des/camcat/Data/Catalogs/Others/jma_cat_2010_2013_update20130329_sel.dat";
	//char cat_file[]="/home/des/camcat/Data/Catalogs/Others/jma_cat_2011_update20130329.dat";
	char fname[120];
	//char fore_file[]="input/other/tohoku_bgrate_nonuniform_3to39.dat";
	char fore_file[]="input/other/tohoku_template_uniform.dat";
	//char fore_file[]="input/other/darf_temp_0.05.dat";
	struct crust crst;
	struct tm reftime;
	//double t0=-434.614525;	//from 01/01/2010 (start of catalog).
	//double t0=-79.0;	//from 01/01/2011 (start of catalog).
	double Mmain=8.5;
	double dR=50, dZ=50;
	int ord=2, err=0;
	double smoothing=5.0;	//min. distance used for smoothing.
	double Mc;
	double *rategrid=NULL;


	//sscanf("2012-09-03T16:35:42Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	sscanf("2012-03-11T14:46:18Z", "%d-%d-%dT%d:%d:%dZ", &(reftime.tm_year), &(reftime.tm_mon), &(reftime.tm_mday), &(reftime.tm_hour), &(reftime.tm_min), &(reftime.tm_sec));
	reftime.tm_year-=1900;
	reftime.tm_mon-=1;
	reftime.tm_isdst=0;

	err=read_crust(fore_file, NULL, &crst, 100.0, 100.0, 0);
	if (err) return;

	//sprintf(fname,"test/bgrate_toho2010_%.2f.dat",Mcuts[i]);
	//sprintf(fname,"test/bgrate_darf.dat");
	sprintf(fname,"test/bgrate_tohoB.dat");
	err=background_rate(cat_file, &crst, reftime, Mmain, &Mc, NULL, &rategrid, dR, dZ, smoothing, ord);
	crst.nmags=1;
	crst.GRmags=dvector(1,1);
	crst.GRmags[1]=1.0;
	csep_forecast(fname, crst, rategrid, 0);

	return;
}

