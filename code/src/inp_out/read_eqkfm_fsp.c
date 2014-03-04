/*
 * read_eqkfm.c
 *
 *  Created on: Jan 16, 2013
 *      Author: camcat
 */

#include "read_eqkfm.h"

int read_fsp_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF_out){

	int ff;
	int nf;
	FILE *fin;
	double lat0, lon0, dep0, str0, dip0, rake0;	//epicenter.
	double l0, d0, Mw, value;
	char errmsg[300];
	int ndi_nst, err=0;
	double dL, dW;

	sprintf(errmsg, "Could not read input slip model %s (read_fsp_eqkfm).", fname);

	fin = fopen(fname, "r");
	if (fin==NULL){
		if (verbose_level>0) printf("Invalid input file (%s).\n", fname);
		if (flog) fprintf(flog,"Invalid input file passed to read_pscmp_eqkfm (%s).\n", fname);
		return (1);
	}

	if (!(ff=next_separator(fin,"FINITE-SOURCE RUPTURE MODEL"))) error_quit(errmsg);
	find_key(fin, "LAT", &lat0);
	find_key(fin, "LON", &lon0);
	find_key(fin, "DEP", &dep0);
	find_key(fin, "LEN", &l0);
	find_key(fin, "WID", &d0);
	find_key(fin, "STRK", &str0);
	find_key(fin, "DIP", &dip0);
	find_key(fin, "RAKE", &rake0);
	find_key(fin, "Mw", &Mw);
	if (!(ff=next_separator(fin,"inversion-related parameters"))) error_quit(errmsg);
	find_key(fin, "Nsg", &value);
	nf= (int) value;

	if (NF_out) *NF_out=nf;
	if (eqfm_out){

		*eqfm_out=eqkfm_array(0,nf-1);

		if (nf==1) {
			find_key(fin, "Nx", &value);
			(*eqfm_out)[0].np_st=(int) value;
			find_key(fin, "Nz", &value);
			(*eqfm_out)[0].np_di=(int) value;
			find_key(fin, "Dx", &dL);
			find_key(fin, "Dz", &dW);
			for (int i=1; i<=3; i++) next_separator(fin,"");
			find_key(fin, "Nsbfs", &value);
			ndi_nst=(int) value;
//			next_separator(fin,"");
			if ((*eqfm_out)[0].np_st*(*eqfm_out)[0].np_di!=ndi_nst){
				if (verbose_level) printf("Error: geometry of slip model %s not understood. Exit.\n", fname);
				if (flog) fprintf(flog, "Error: geometry of slip model %s not understood. Exit.\n", fname);
				fclose(fin);
				return 1;
			}
			(*eqfm_out)[0].whichfm=1;
			(*eqfm_out)[0].str1=str0;
			(*eqfm_out)[0].dip1=dip0;
			(*eqfm_out)[0].rake1=rake0;
			(*eqfm_out)[0].L=l0;
			(*eqfm_out)[0].W=d0;
			(*eqfm_out)[0].depth=dep0;
			(*eqfm_out)[0].lat=lat0;
			(*eqfm_out)[0].lon=lon0;

			(*eqfm_out)[0].pos_d=dvector(1,ndi_nst);
			(*eqfm_out)[0].pos_s=dvector(1,ndi_nst);
			(*eqfm_out)[0].slip_str=dvector(1,ndi_nst);
			(*eqfm_out)[0].slip_dip=dvector(1,ndi_nst);

			err+=read_slipvalues(fin, *eqfm_out);

		}


		else {
			next_separator(fin,"MULTISEGMENT MODEL");
			next_separator(fin,"");
			for (int f=0; f<nf; f++){
				(*eqfm_out)[f].whichfm=1;
				find_key(fin, "STRIKE", &(*eqfm_out)[f].str1);
				find_key(fin, "DIP", &(*eqfm_out)[f].dip1);
				find_key(fin, "RAKE", &(*eqfm_out)[f].rake1);
				find_key(fin, "LEN", &(*eqfm_out)[f].L);
				find_key(fin, "WID", &(*eqfm_out)[f].W);
				find_key(fin, "Z2top", &(*eqfm_out)[f].depth);
				find_key(fin, "LAT", &(*eqfm_out)[f].lat);
				find_key(fin, "LON", &(*eqfm_out)[f].lon);
				find_key(fin, "Nsbfs", &value);
				ndi_nst=(int) value;
				find_key(fin, "Dx", &dL);
				find_key(fin, "Dz", &dW);

				//np_di, np_st are also given in the file (Nx, Nz), but only once -> what about models in which each subfault has different no of patches?
				(*eqfm_out)[f].np_di=(int) ((*eqfm_out)[f].W/dW);
				(*eqfm_out)[f].np_st=(int) ((*eqfm_out)[f].L/dL);
				(*eqfm_out)[f].pos_d=dvector(1,ndi_nst);
				(*eqfm_out)[f].pos_s=dvector(1,ndi_nst);
				(*eqfm_out)[f].slip_str=dvector(1,ndi_nst);
				(*eqfm_out)[f].slip_dip=dvector(1,ndi_nst);

				if ((*eqfm_out)[f].np_st*(*eqfm_out)[f].np_di!=ndi_nst){
					if (verbose_level) printf("Error: geometry of slip model %s not understood. Exit.\n", fname);
					if (flog) fprintf(flog, "Error: geometry of slip model %s not understood. Exit.\n", fname);
					return 1;
				}
				err+=read_slipvalues(fin, (*eqfm_out)+f);
			}
		}
	}

	fclose(fin);
	return (err!=0);

}

void track_position(long *pos, int NP, FILE* fin){	//todo move to more general file.

	for (int i=NP; i>0; i--) pos[i]=pos[i-1];
	pos[0]=ftell(fin);
	return;
}

int next_separator(FILE *fin, char *string){

	char comment1[]="% -";
	char comment2[]="%--";
	int Nchar_long=500;
	int nchar=strlen(string), off;
	int found=0, is_sep;
	char line[Nchar_long];
	while (!feof(fin) && !found){
		is_sep=0;
		while(!feof(fin) && !is_sep){
			line[0]=0;
			while(line[0]!=comment1[0]) fgets(line,Nchar_long,fin);
			if (!strncmp(line,comment1,3) || !strncmp(line,comment2,3)) is_sep=1;
		}
		off=0;
		while(off+nchar<strlen(line)){
			if (!strncmp(line+off,string,nchar)) {
				found=1;
				break;
			}
			off++;
		}
	}

	return found;
}

int find_key(FILE *fin, char *string, double *value){

	char comment1[]="% -";
	char comment2[]="%--";
	int Nchar_long=500;
	int nchar=strlen(string), off;
	int found=0;
	long int pos0=ftell(fin);
	char line[Nchar_long];

	strncpy(line, comment1, 3);
	while (!feof(fin) && (!strncmp(line,comment1,3) || !strncmp(line,comment2,3) || line[0]!=comment1[0])) fgets(line,Nchar_long,fin);

	while (!feof(fin) && !found){
		if (!strncmp(line,comment1,3) || !strncmp(line,comment2,3) || line[0]!=comment1[0]) break;
		off=0;
		while(off+nchar<strlen(line)){
			if (!strncmp(line+off,string,nchar)) {
				found=1;
				sscanf(line+off,"%*s = %lf", value);
				break;
			}
			off++;
		}
		fgets(line,Nchar_long,fin);
	}

	fseek (fin, pos0, SEEK_SET);
	return found;
}

int scan_nth(char *string, int n, double *result){

	int numread=1;
	char line[1000];
	char linedum[1000];

	sprintf(line,"%s", string);

	for (int i=1; i<=n; i++){

		sprintf(linedum,"");
		numread=sscanf(line,"%lf %[^\n]", result, linedum);
		if (numread<1) return 1;
		else sprintf(line,"%s",linedum);
	}

	return 0;
}

int read_slipvalues(FILE *fin, struct eqkfm *eqfm){
	//stops reading at next comment ("%").

	char comment[]="%", comm=comment[0];
	int Nchar_long=1500;
	char line[Nchar_long], linedum[Nchar_long], key[10];
	double lat, lon, dep, slip, rake, lattop, lontop;
	double dipVx, dipVy, dx, dy;
	int i=0, ntimes=0;
	long pos[4];
	double strike_toll=15*DEG2RAD;
	int nr_slip=0, nr_rake=0;
	int rakefound=0, slipfound=0, numread=1;

	double 	lat0=(*eqfm).lat,		//
			lon0=(*eqfm).lon,		//
			dep0=(*eqfm).depth,		//
			strike=(*eqfm).str1,	//
			dip=(*eqfm).dip1,		//
			rake0=(*eqfm).rake1,		//
			dw=(*eqfm).W/(*eqfm).np_di;
	double 	*pos_str=(*eqfm).pos_s,	//
			*pos_dip=(*eqfm).pos_d,	//
			*slip_str=(*eqfm).slip_str,	//
			*slip_dip=(*eqfm).slip_dip;

	strike*=DEG2RAD;
	dip*=DEG2RAD;

	dipVx=cos(dip)*sin(strike+pi/2);
	dipVy=cos(dip)*cos(strike+pi/2);

	line[0]=comm;
	while (line[0]==comm) {
		fgets(line,Nchar_long,fin);
		track_position(pos, 3,fin);
		ntimes+=1;
	}

	if (ntimes<3){
		if (verbose_level) printf("**Warning: could not find columns labels in input file, will assume the following: LAT LON X Y Z SLIP RAKE (read_slipvalues)**\n");
		if (flog) fprintf(flog, "**Warning: could not find columns labels in input file, will assume the following: LAT LON X Y Z SLIP RAKE (read_slipvalues)**\n");
	}
	fseek (fin , pos[3] , SEEK_SET);
	fgets(line,Nchar_long,fin);

	while (numread>0 && (!rakefound || !slipfound)){
		sprintf(linedum,"");
		numread=sscanf(line,"%s %[^\n]", key, linedum);
		sprintf(line,"%s",linedum);
		if (!strncmp(key, "SLIP", 4)) slipfound=1;
		else if (!slipfound) nr_slip+=1;
		if (!strncmp(key, "RAKE", 4)) rakefound=1;
		else if (!rakefound) nr_rake+=1;
	}

	if (!rakefound){
		if (verbose_level>1) printf("** Warning: could not find field named RAKE. Will use global value (%.3lf). **\n", rake0);
		if (flog) {
			fprintf(flog, "** Warning: could not find field named RAKE. Will use global value (%.3lf). **\n", rake0);
			fflush(flog);
		}
	}
	if (!slipfound){
		if (verbose_level) printf("** Error: could not find field named SLIP. Exiting. **\n");
		if (flog) {
			fprintf(flog, "** Error: could not find field named SLIP. Exiting. **\n");
			fflush(flog);
		}
		return (1);
	}

	fgets(line,Nchar_long,fin);
	fgets(line,Nchar_long,fin);

	while (!feof(fin)){
		i+=1;
		sscanf(line, "%lf %lf %*f %*f %lf", &lat, &lon, &dep);
		scan_nth(line, nr_slip, &slip);
		if (rakefound) scan_nth(line, nr_rake, &rake);
		else rake=rake0;
		pos_dip[i]=(dep-dep0)/sin(dip)+0.5*(dw);
		lontop=lon-RAD2DEG*(pos_dip[i]-0.5*(dw))*dipVx/(Re*cos(DEG2RAD*lat));
		lattop=lat-RAD2DEG*(pos_dip[i]-0.5*(dw))*dipVy/Re;
		dy=Re*(lattop-lat0)*DEG2RAD;
		dx=(Re*cos(DEG2RAD*lat))*(lontop-lon0)*DEG2RAD;
		if ((dy>0.1) && (fabs(fmod(atan(dx/dy)+pi,pi)-fmod(strike,pi))>strike_toll)){
			if (verbose_level>1){
				printf(" **Warning: Epicenter is not on the fault! Will use a new epicenter ** \n");
				printf("\tOld epicenter: [lat, lon]= [%.3lf, %.3lf]\n", lat0,lon0);
			}
			if (flog){
				fprintf(flog," **Warning: Epicenter is not on the fault! Will use a new epicenter ** \n");
				fprintf(flog,"\tOld epicenter: [lat, lon]= [%.3lf, %.3lf]\n", lat0,lon0);
			}
			lat0=(*eqfm).lat=lattop;
			lon0=(*eqfm).lon=lontop;
			if (verbose_level>1) printf("\tNew epicenter: [lat, lon]= [%.3lf, %.3lf]\n", lat0,lon0);
			if (flog){
				fprintf(flog, "\tNew epicenter: [lat, lon]= [%.3lf, %.3lf]\n", lat0,lon0);
				fflush(flog);
			}
			fseek (fin , pos[3] , SEEK_SET);
			for (int h=1; h<=3; h++) fgets(line,Nchar_long,fin);
			i=0;
			continue;
		}
		pos_str[i]=sqrt(dx*dx+dy*dy)*sign(dx/sin(strike));
		slip_str[i]=slip*cos(DEG2RAD*rake);
		slip_dip[i]=-slip*sin(DEG2RAD*rake);
		fgets(line,Nchar_long,fin);
		if (line[0]==comm) break;
	}
	if (i>(*eqfm).np_di*(*eqfm).np_st) {
		if (verbose_level) printf("**Error: wrong number of lines in slip model file, or string 'line' is too short (read_slipvalues). Exiting. \n**");
		if (flog) fprintf(flog, "**Error: wrong number of lines in slip model file, or string 'line' is too short (read_slipvalues). Exiting. \n**");
		return 1;
	}

	return 0;
}

