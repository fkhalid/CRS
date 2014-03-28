/*
 * readfocmec.c

 *
 *  Created on: Feb 17, 2012
 *      Author: camcat
 */


#include "read_focmec.h"
#include "mpi.h"

int readmultiplefocmec(char **focmecfiles, int nofiles, char *which_format,
					  struct crust crst, double border, double dz, double dDCFS,
					  struct tm reftime, double t0, double t1, double tfocmec,
					  double mag, double ***focmec,	int **firstelements, int *NFM,
					  int *NFM_timesel, struct eqkfm **eqkfm,int sel, int fm2) {
// *firstelements should not be initialized (is done inside this function).

	double **focmectemp;
	int nfmtemp, nfmtemp2, cl;
	int ntotmax=0, err=0, nfm_sofar=0, nfm_sofar2=0;
	struct eqkfm *eqkfmtemp;

	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	if(procId == 0) {
		for(int n=0; n<nofiles; n++) {
			cl=countline(focmecfiles[n]);

			if (cl>=0) {
				ntotmax+=countline(focmecfiles[n]);
			}
			else {
				fileError = 1;
				break;
			}
		}
	}
	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif
	if(fileError) {
		return 1;
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&ntotmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (fm2) ntotmax*=2;

	if(focmec) {
		*focmec = dmatrix(1,4,1,ntotmax);
	}
	if (firstelements) *firstelements=ivector(0,nofiles);
	if (eqkfm) *eqkfm=eqkfm_array(0,ntotmax-1);

	for (int n=0; n<nofiles; n++) {
		err += readfocmec(focmecfiles[n], which_format, crst, border, dz, dDCFS, reftime, t0,
						  t1, tfocmec, mag, (focmec)? &focmectemp : NULL, &nfmtemp, &nfmtemp2,
						  (eqkfm)? &eqkfmtemp : NULL, sel, fm2);

		if (focmec && !nfmtemp) {
			if(procId == 0) {
				if (verbose_level) printf("**Warning: no focal mechanisms selected from file %s. (readmultiplefocmec).**\n", focmecfiles[n]);
				if (flog) {
					fprintf(flog, "**Warning: no focal mechanisms selected from file %s. (readmultiplefocmec).**\n", focmecfiles[n]);
					fflush(flog);
				}
			}
		}

		if (firstelements) (*firstelements)[n]=nfm_sofar+1;
		if (focmec) {
			for (int ns=1; ns<=4; ns++){
				for (int n=1; n<=nfmtemp; n++) (*focmec)[ns][n+nfm_sofar]=focmectemp[ns][n];
			}
			nfm_sofar+=nfmtemp;
		}
		if (eqkfm){
			for (int i=0; i<nfm_sofar2+nfmtemp2; i++) copy_eqkfm_all(eqkfmtemp[i], (*eqkfm)+nfm_sofar2+i);
			nfm_sofar2+=nfmtemp2;
		}
	}

	if (firstelements) (*firstelements)[nofiles]=nfm_sofar+1;
	if (NFM) *NFM=nfm_sofar;
	if (NFM_timesel) *NFM_timesel=nfm_sofar2;

	return (err>0);
}

int readfocmec(char *focmecfile, char *which_format, struct crust crst,
			   double border, double dz, double dDCFS, struct tm reftime,
			   double t0, double t1, double tfocmec, double mag,
			   double ***focmec, int *NFM, int *NFM_timesel,
			   struct eqkfm **eqkfm,int sel, int fm2) {
// border, dz indicate extra volume to be considered for spatial selection.
// magnitude selection only applies to sources, not to contents of focmec.
// sel is flag to specify if spatial selection should be done.
// fm2 is flag to specify if both mechanisms should be used.
// focmec will contain all focal mechanisms in the relevant area, and time t<=tfocmec.  size [1...NC,1...NFM]
// eqkfm will only contain those between time t0, t1 (i.e. to be used as sources). size [0...NFM_timesel]
// if focmec=NULL or eqkfm=NULL, they will be ignored.

	// [Fahad] Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	//check if file exists and can be opened.
	if(procId == 0) {
		if (flog) {
			fprintf(flog, "\nReading focal mechanisms from file %s (format: %s).\n", focmecfile, which_format);
			if (fm2) fprintf(flog, "Using both focal mechanisms.\n");
			else fprintf(flog, "Using only first focal mechanism.\n");
			fflush (flog);
		}
	}

	FILE *fin;

	if(procId == 0) {
		fin = fopen(focmecfile,"r");
		if(fin == NULL) {
			if (verbose_level>1) printf("** Error: could not open focal mechanisms catalog %s (readfocmec.c). **\n", focmecfile);
			if (flog) {
				fprintf(flog, "Error: could not open focal mechanisms catalog (readfocmec.c).\n");
				fflush (flog);
			}

			fileError = 1;
		}
		else {
			fclose(fin);
		}
	}
	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif
	if(fileError) {
		return (1);
	}

	int NFMmax;
	double **focmec0;
	int *selected, *selectedsources, p;
	int NC;
	int NFM2=0, NFM2sources=0;
	long NFM0;
	int ignore_time=0;
	int err=0;
	int H, lat_col, lon_col, dep_col, mag_col, str_col, dip_col, rake_col, time_col;	//H=no. of header lines.
	struct tm ev;
	char timestr[20];
	double slip;
	double *times;
	double 	latmin=crst.latmin-180*border/(Re*PI), \
			latmax=crst.latmax+180*border/(Re*PI), \
			lonmin=crst.lonmin-180*border/(Re*PI*cos(crst.lat0*PI/180)), \
			lonmax=crst.lonmax+180*border/(Re*PI*cos(crst.lat0*PI/180)), \
			depthmin=fmin(0.0, crst.depmin-dz), \
			depthmax=crst.depmax+dz;

	int *taper_all;
	taper_all=ivector(1,4);
	for (int i=1; i<=4; i++) taper_all[i]=1;

	if(procId == 0) {
		NFMmax = (fm2==1)? 2*countline(focmecfile) : countline(focmecfile);
		NC = countcol(focmecfile);
	}
	#ifdef _CRS_MPI
		MPI_Bcast(&NFMmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NC, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	// FIXME: [Fahad] Any 'return' statement other than the one at the end of the function
	//		   should be prominent enough so that it is not missed by another
	//		   person looking at the code.
//	if (NFMmax<0) return 1;
	if(NFMmax < 0) {
		return 1;
	}

	if(procId == 0) {
		if (verbose_level>1) printf("Reading catalog of focal mechanisms...");
	}

	if (strcmp(which_format,"7col")==0){
		lat_col=1;
		lon_col=2;
		dep_col=3;
		mag_col=4;
		str_col=5;
		dip_col=5;
		rake_col=6;
		time_col=0;
		H=0;

		if (eqkfm) {
			if(procId == 0) {
				if (verbose_level>0) printf("** Warning: non null value of eqkfm passed to readfocmec, but structure can not be set up since no time information contained in input file.**\n");
			}
			*eqkfm=NULL;
			ignore_time=1;
		}
	}

	else{
		if (strcmp(which_format,"CSEP")==0){
			lat_col=3;
			lon_col=4;
			dep_col=14;
			mag_col=12;
			str_col=5;	//str, dip, rake have also columns 8, 9, 10.
			dip_col=6;
			rake_col=7;
			time_col=2;
			H=1;
		}
		else {
			if(procId == 0) {
				if (verbose_level>0) printf ("**Error: format identifier not recognized (readfocmec.c). Exiting.** \n");
				if (flog) fprintf(flog,"Error: format identifier not recognized (readfocmec.c). Exiting.\n");
			}

			return (1);
		}
	}

	focmec0=dmatrix(1,NC,1,NFMmax);
	selected=ivector(1,NFMmax);
	selectedsources=ivector(1,NFMmax);
	times=dvector(1,NFMmax);

	if(procId == 0) {
		err = read_matrix(focmecfile, NC, H, focmec0, &NFM0);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NFM0, 1, MPI_LONG, 0, MPI_COMM_WORLD);

		long nrl=1, nrh=NC, ncl=1, nch=NFMmax;
		long nrow=nrh-nrl+1, ncol=nch-ncl+1;
		MPI_Bcast(focmec0[nrl], nrow*ncol+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	if(procId == 0) {
		for (int p=1; p<=NFM0; p++){
			if ((sel==0) | ((focmec0[lat_col][p]>=latmin && focmec0[lat_col][p]<=latmax && focmec0[lon_col][p]>=lonmin && focmec0[lon_col][p]<=lonmax && focmec0[dep_col][p]>=(depthmin-1) && focmec0[dep_col][p]<=(depthmax+1)))){
				if (!ignore_time){
					sprintf(timestr, "%.0lf", focmec0[time_col][p]);
					sscanf(timestr, "%4d%2d%2d%2d%2d%2d", &(ev.tm_year), &(ev.tm_mon), &(ev.tm_mday), &(ev.tm_hour), &(ev.tm_min), &(ev.tm_sec));
					ev.tm_year-=1900;
					ev.tm_mon-=1;
					times[p]=difftime(mktime(&ev),mktime(&reftime))*SEC2DAY;
				}		

				if (focmec){
					if (ignore_time || (times[p]<=tfocmec)){
						NFM2+=1;
						selected[NFM2]=p;
					}
				}

				if (eqkfm){
					if (ignore_time || (times[p]>=t0 && times[p]<=t1 && focmec0[mag_col][p]>=mag)){
						NFM2sources+=1;
						selectedsources[NFM2sources]=p;
					}
				}
			}
		}

		if (flog) {
			fprintf(flog, "%d events selected as sample of focal planes, %d selected as sources.\n", NFM2, NFM2sources);
			fflush(flog);
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NFM2, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NFM2sources, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(times, NFMmax+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(selected, NFMmax+1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(selectedsources, NFMmax+1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	//-------------fill in matrix of focal mechanisms-----------------//

	if (NFM) *NFM= (focmec)? (fm2 ? 2*NFM2 : NFM2) : 0;
	if (focmec){
		*focmec=dmatrix(1,4,1,*NFM);
		for (int p0=1; p0<=NFM2; p0++){
			p=selected[p0];
			(*focmec)[1][p0]=focmec0[str_col][p];
			(*focmec)[2][p0]=focmec0[dip_col][p];
			(*focmec)[3][p0]=focmec0[rake_col][p];
			(*focmec)[4][p0]= (ignore_time) ? -1e30 : times[p];

			if (fm2==1){
				(*focmec)[1][p0+NFM2]=focmec0[str_col+3][p];
				(*focmec)[2][p0+NFM2]=focmec0[dip_col+3][p];
				(*focmec)[3][p0+NFM2]=focmec0[rake_col+3][p];
				(*focmec)[4][p0+NFM2]=(ignore_time) ? -1e30 : times[p];
			}
		}
	}

	//-------------fill in catalog of sources with foc mec---------------//

	if (NFM_timesel) *NFM_timesel= (eqkfm)? NFM2sources : 0;
	if (eqkfm){
		*eqkfm=eqkfm_array(0,NFM2sources-1);
		for (int p0=0; p0<NFM2sources; p0++){		//todo parallel
			p=selectedsources[p0+1];
			(*eqkfm)[p0].t=times[p];
			(*eqkfm)[p0].lat=focmec0[lat_col][p];
			(*eqkfm)[p0].lon=focmec0[lon_col][p];
			(*eqkfm)[p0].depth=focmec0[dep_col][p];
			(*eqkfm)[p0].mag=focmec0[mag_col][p];
			(*eqkfm)[p0].index_cat= 0;
			(*eqkfm)[p0].str1=focmec0[str_col][p];
			(*eqkfm)[p0].dip1=focmec0[dip_col][p];
			(*eqkfm)[p0].rake1=focmec0[rake_col][p];
			if ((*eqkfm)[p0].rake1<0) (*eqkfm)[p0].rake1+=360;
			if (fm2){
				(*eqkfm)[p0].str2=focmec0[str_col][p+3];
				(*eqkfm)[p0].dip2=focmec0[dip_col][p+3];
				(*eqkfm)[p0].rake2=focmec0[rake_col][p+3];
				if ((*eqkfm)[p0].rake2<0) (*eqkfm)[p0].rake2+=360;
				(*eqkfm)[p0].whichfm=0;
			}
			else (*eqkfm)[p0].whichfm=1;

			(*eqkfm)[p0].is_slipmodel=1;
			(*eqkfm)[p0].np_st=1;
			(*eqkfm)[p0].np_di=1;
			(*eqkfm)[p0].pos_s=dvector(1,1);	//location of patches within fault; [0], [0] for single patch events.
			(*eqkfm)[p0].pos_d=dvector(1,1);
			(*eqkfm)[p0].taper=taper_all;
			(*eqkfm)[p0].pos_s[1]=0;	//location of patches within fault; [0], [0] for single patch events.
			(*eqkfm)[p0].pos_d[1]=0;
			if ((*eqkfm)[p0].whichfm){
				(*eqkfm)[p0].slip_str=dvector(1,1);
				(*eqkfm)[p0].slip_dip=dvector(1,1);
			}
			else {
				(*eqkfm)[p0].slip_str=dvector(1,2);
				(*eqkfm)[p0].slip_dip=dvector(1,2);
			}
			err+=find_gridpoints_d(crst.y, crst.x, crst.depth, (int *) 0, 0, crst.N_allP, (*eqkfm)[p0].y, (*eqkfm)[p0].x, (*eqkfm)[p0].depth,  (*eqkfm)[p0].mag, dDCFS,  &((*eqkfm)[p0].nsel), &((*eqkfm)[p0].selpoints));
			WellsCoppersmith((*eqkfm)[p0].mag, (*eqkfm)[p0].rake1, &((*eqkfm)[p0].L), &((*eqkfm)[p0].W), &slip);
			slip=(*eqkfm)[p0].tot_slip=pow(10,(1.5*((*eqkfm)[p0].mag+6)))*(1.0/(crst.mu*pow(10,12)*(*eqkfm)[p0].W*(*eqkfm)[p0].L));
			if ((*eqkfm)[p0].depth<0.5*(*eqkfm)[p0].W*sin(DEG2RAD*(*eqkfm)[p0].dip1)) (*eqkfm)[p0].depth=0.5*(*eqkfm)[p0].W*sin(DEG2RAD*(*eqkfm)[p0].dip1);
			(*eqkfm)[p0].slip_str[1]=slip*cos(DEG2RAD*(*eqkfm)[p0].rake1);
			(*eqkfm)[p0].slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfm)[p0].rake1);

			//by convention, slip_xxx[2] contains the slip for second foc mech (for single patch events only!).
			if (!(*eqkfm)[p0].whichfm){
				if ((*eqkfm)[p0].depth<0.5*(*eqkfm)[p0].W*sin(DEG2RAD*(*eqkfm)[p0].dip2)) (*eqkfm)[p0].depth=0.5*(*eqkfm)[p0].W*sin(DEG2RAD*(*eqkfm)[p0].dip2);
				(*eqkfm)[p0].slip_str[2]=slip*cos(DEG2RAD*(*eqkfm)[p0].rake2);
				(*eqkfm)[p0].slip_dip[2]=-slip*sin(DEG2RAD*(*eqkfm)[p0].rake2);
			}
	    	latlon2localcartesian((*eqkfm)[p0].lat, (*eqkfm)[p0].lon, crst.lat0, crst.lon0, &((*eqkfm)[p0].y), &((*eqkfm)[p0].x));

		}
	}

	free_dmatrix(focmec0,1,NC,1,NFMmax);
	free_ivector(selected,1,NFMmax);
	free_ivector(selectedsources,1,NFMmax);
	free_dvector(times,1,NFMmax);

	if(procId == 0) {
		if (verbose_level>1) printf("done.\n");
	}

	return (err!=0);

}

void select_fm_time(double **focmec, int *NFM, double Tstart){

	int ntot=0;

	for (int n=1; n<=*NFM; n++){
		if (focmec[4][n]<Tstart) {
			ntot+=1;
			for (int i=1; i<=4; i++) focmec[i][ntot]=focmec[i][n];
		}
	}

	*NFM=ntot;
}
