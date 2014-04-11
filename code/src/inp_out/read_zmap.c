/*
 * main.c
 *
 *  Created on: Sep 15, 2013
 *      Author: camilla
 */

#include "read_zmap.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int readZMAP (struct catalog *cat, struct eqkfm **eqfm, int *Ntot, char *file,
			  struct crust crst, struct tm reftime, double t0s, double t1s,
			  double t0c, double t1c, double Mmain, double tw, double border,
			  double extra_d, double dDCFS, int findgridpoints) {
/* t0s,t1s: for considering sources.
 * t0c,t1c: for catalog (LL periods).
 * a time window of length tw will be discarded after each event of magnitude >=Mmain in calculating Mc, b values. (these events should not be used in the LL calculations).
 * all are given in days from reftime.
 * border (in km): extra distance to be considered for sources.
 * cat, eqfm can be given as NULL, and will be skipped.
 * findgridpoints = flag indicating if grid points corresponding to events should be calculated (expensive).
 * Otherwise, everything in these structures is initialized.
 *
 *  * */

	// [Fahad] Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double std_merr=0.1, std_verr=5.0, std_herr=4.0;	//todo read from somewhere!
	int line_length=2000;
	int hh, lines=0, valid=0, empty=0, missing_values=0;
	int * old2new;
	FILE *fin;
	char line[line_length], *st;
	int no_expected_columns=13;	//ZMAP format as used by CSEP?
	int Z;
	double Mc_offset=0.3, dM=0.9;	//todo don't hardwide Mc_offset. todo dM should be passed.
	int cut_sd=3.0;	//no. of s.dev. for cutting off gaussian.	//todo read from somewhere?
	double t_last_large;
	int k;
	struct tm ev;
	int lon_out_of_range, lat_out_of_range, date_out_of_range, time_out_of_range, mag_out_of_range, dep_out_of_range;

	int mon, day, hour, min, sec;
	double fyear, fmon, fday, fhour, flmin, fsec;
	double year, this_year;
	time_t t;

	int eq1=0, eq2=0, eq;
	int errP=0;
	int *seleq1, *seleq2, *catindex;	//catindex: indices of events in catalog, to be copied into eqkfm.
	double SD, SDd, SDlat, SDlon, f=1.0;
	double lat0l, lat1l, lon0l, lon1l, dep0l, dep1l;
	double x,y;
	double Mc;
	double 	*xgrid=crst.x, \
			*ygrid=crst.y, \
			*depgrid=crst.depth, \
			*dAgrid=crst.dAgrid;
	int N=crst.N_allP;

	if(procId == 0) {
		Z = countline(file)+1;
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&Z, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	// FIXME: [Fahad] Any 'return' statement other than the one at the end of the function
	//		   should be prominent enough so that it is not missed by another
	//		   person looking at the code.
//	if (Z<=0) return 1;
	if (Z <= 0) {
		return 1;
	}

	t=time(NULL);
	this_year=(*localtime(&t)).tm_year+1900;

	double 	lon0=crst.lonmin, \
			lon1=crst.lonmax, \
			lat0=crst.latmin, \
			lat1=crst.latmax, \
			dep0=crst.depmin, \
			dep1=crst.depmax;

	double 	*lat=dvector(1,Z), \
			*lon=dvector(1,Z), \
			*mag=dvector(1,Z), \
			*mag2=dvector(1,Z), \
			*dep=dvector(1,Z), \
			*herr=dvector(1,Z), \
			*verr=dvector(1,Z), \
			*merr=dvector(1,Z), \
			*times=dvector(1,Z);  //days from start of catalog (t1c)?

	setenv("TZ", "UTC", 1);
	
	if(procId == 0) {	
		fin = fopen(file,"r");

		if(fin == NULL) {
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		if(procId == 0) {
			if (verbose_level) printf("**Error: unable to open input file %s (readZMAP).**\n", file);
			if (flog) fprintf(flog,"**Error: unable to open input file %s (readZMAP).**\n", file);
		}

		return (1);
	}


	//------------------------------read entire catalog:-------------------------//

	int exp_not=0;
	char ex[]="e";

	if(procId == 0) {
		while (!feof(fin)) {
			lines+=1;
			st=fgets(line,line_length,fin);

			if (lines==1) {
				for (int i=0; i<strlen(line) && !exp_not; i++) {
					if (line[i]==ex[0]) exp_not=1;
				}
			}

			//initialize to out of range values (so will notice if columns are missing):
			lon[valid+1]=999;
			lat[valid+1]=999;
			year=0;
			mon=12;
			day=32;
			mag[valid+1]=999;
			dep[valid+1]=-999;
			hour=25;
			min=62;
			sec=62.0;

			//initialize to standard values (in case they are not given in catalog).
			merr[valid+1]=std_merr;
			verr[valid+1]=std_verr;
			herr[valid+1]=std_herr;

			if (exp_not) {
				hh=sscanf(line, "%25le %25le %25le %25le %25le %25le %25le %25le %25le %25le %25le %25le %25le",
							lon+valid+1, lat+valid+1, &fyear, &fmon, &fday, mag+valid+1, dep+valid+1, &fhour, &flmin, &fsec, herr+valid+1, verr+valid+1, merr+valid+1);

				year=(int) fyear;
				mon=(int) fmon;
				day=(int) fday;
				hour=(int) fhour;
				min= (int) flmin;
				sec= (int) fsec;
			}
			else {
				hh=sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",	lon+valid+1, lat+valid+1, &fyear, &fmon, &fday, mag+valid+1, dep+valid+1, &fhour, &flmin, &fsec, herr+valid+1, verr+valid+1, merr+valid+1);
				year=(int) fyear;
				mon=(int) fmon;
				day=(int) fday;
				hour=(int) fhour;
				min= (int) flmin;
				sec= (int) fsec;
			}

			if (fabs(herr[valid+1])<tol0) herr[valid+1]=std_herr;
			if (fabs(verr[valid+1])<tol0) verr[valid+1]=std_herr;
			if (fabs(merr[valid+1])<tol0) merr[valid+1]=std_herr;

			ev.tm_year=floor(year)-1900;
			ev.tm_mon=mon-1;
			ev.tm_mday=day;
			ev.tm_isdst=0;
			ev.tm_hour=hour;
			ev.tm_min=min;
			ev.tm_sec=(int) sec;

			times[valid+1]=difftime(mktime(&ev),mktime(&reftime))*SEC2DAY;

			if (st==NULL) empty+=1;
			else {
				lon_out_of_range = lat_out_of_range = date_out_of_range = time_out_of_range = mag_out_of_range = dep_out_of_range = 0;

				if (lon[valid+1]<-180 || lon[valid+1]>360) lon_out_of_range=1;
				if (lat[valid+1]<-90 || lat[valid+1]>90) 	 lat_out_of_range=1;
				if (year<1000 || year>this_year ||	mon<1 || mon>12 || day<1 || day>31) date_out_of_range=1;
				if (hour<0 || hour>24 ||	min<0 || min>59 || sec<0 || sec>61) time_out_of_range=1;
				if (mag[valid+1]<-10 || mag[valid+1]>12) mag_out_of_range=1;
				if (dep[valid+1]<0 || dep[valid+1]>2000) dep_out_of_range=1;

				if (lon_out_of_range || lat_out_of_range || date_out_of_range || time_out_of_range || mag_out_of_range || dep_out_of_range) {
					missing_values+=1;
					if (hh!=no_expected_columns || missing_values){
						if (verbose_level>1) {
							printf("** Warning: line %d has following columns out of range:", lines);
							if (lon_out_of_range) printf("lon, ");
							if (lat_out_of_range) printf("lat, ");
							if (dep_out_of_range) printf("dep, ");
							if (mag_out_of_range) printf("mag, ");
							if (date_out_of_range) printf("date, ");
							if (time_out_of_range) printf("time, ");
							printf(" and will be skipped. ** \n");
						}
						if (flog) {
							fprintf(flog, "Warning: line %d has following columns out of range:", lines);
							if (lon_out_of_range) fprintf(flog, "lon (%.2lf), ", lon[valid+1]);
							if (lat_out_of_range) fprintf(flog, "lat (%.2lf), ", lat[valid+1]);
							if (dep_out_of_range) fprintf(flog, "dep (%.2lf), ", dep[valid+1]);
							if (mag_out_of_range) fprintf(flog, "mag (%.2lf), ", mag[valid+1]);
							if (date_out_of_range) fprintf(flog, "date (%2d-%2d-%4.lf), ", day, mon, year);
							if (time_out_of_range) fprintf(flog, "time (%2d:%2d:%2d). ", hour, min, sec);
							fprintf(flog, " and will be skipped.\n");
						}
					}
					else valid+=1;
				}
				else valid+=1;
			}
		}
		fclose(fin);

		if (verbose_level > 2) printf("%d of %d lines are valid; %d have missing values; %d are empty\n", valid, lines, missing_values, empty);
		if (flog) fprintf(flog, "%d of %d lines are valid; %d have missing values; %d are empty\n", valid, lines, missing_values, empty);
		else {
			if (valid!=lines && verbose_level>0) {
				printf("** Warning: %d invalid lines skipped in catalog: %s (%d lines have missing values; %d are empty). **\n", missing_values+empty, file, missing_values, empty);
				if (flog) fprintf(flog, "Warning: %d invalid lines skipped in catalog: %s (%d lines have missing values; %d are empty).\n", missing_values+empty, file, missing_values, empty);
			}
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(valid == 0){
		(*cat).Z = 0;

		return(1);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(lat, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(lon, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(mag, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(mag2,  Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(dep, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(herr,  Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(verr,  Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(times, Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	//------------------------------select events:-------------------------//

	//by convention, Mc>20 means that Mc should be calculated y program (later).
	if (!cat || (*cat).Mc>=20) Mc=-10;	//select all events.
	else Mc=(*cat).Mc;		//select events>=(*cat).Mc.

	//define large boundaries (for sources):
	lat0l=lat0-180*border/(Re*PI);
	lat1l=lat1+180*border/(Re*PI);
	lon0l=lon0-180*border/(Re*PI*cos(crst.lat0*PI/180));
	lon1l=lon1+180*border/(Re*PI*cos(crst.lat0*PI/180));
	dep0l=fmin(0.0,dep0-extra_d);
	dep1l=dep1+extra_d;

	//count points in catalog, and sources:
	seleq1=ivector(0,valid);
	catindex=ivector(0,valid);
	seleq2=ivector(0,valid);
	old2new=ivector(0,valid);
	t_last_large=-1e30;
	for (int i=1; i<=valid; i++){
		SD=f*herr[i];
		SDd=f*verr[i];
		SDlat=180*SD/(Re*PI);
		SDlon=180*SD/(Re*PI*cos(0.5*(lat0+lat1)*PI/180));

		if (mag[i]>= Mc-dM && lat[i]+SDlat>=lat0l && lat[i]-SDlat<=lat1l && lon[i]+SDlon>=lon0l && lon[i]-SDlon<=lon1l && dep[i]+SDd>=dep0l && dep[i]-SDd<=dep1l){
			if (times[i]>=t0s && times[i]<=t1s) {
				seleq1[eq1]=i;
				eq1+=1;
			}
			if (mag[i]>= Mc && times[i]>=t0c && times[i]<=t1c && lat[i]+SDlat>=lat0 && lat[i]-SDlat<=lat1 && lon[i]+SDlon>=lon0 && lon[i]-SDlon<=lon1 && dep[i]+SDd>=dep0 && dep[i]-SDd<=dep1) {
				if (times[i]>=t_last_large+tw){
					if (times[i]>=t0s && times[i]<=t1s) catindex[eq1-1]=eq2+1;	//+1 since cat.XX[1...cat.Z].
					seleq2[eq2]=i;
					eq2+=1;
				}
				else if (times[i]>=t0s && times[i]<=t1s) catindex[eq1-1]=0;
			}
			else if (times[i]>=t0s && times[i]<=t1s) catindex[eq1-1]=0;
		}
		if (mag[i]>=Mmain) t_last_large=times[i];
	}

	if (!eq1 && !eq2) {
		if (cat) (*cat).Z=0;
		if (Ntot) *Ntot=0;
		return 1;
	}

	if(procId == 0) {
		if (flog){
			fprintf(flog, "%d events selected for LL inversion. \n", eq2);
			fprintf(flog, "%d events selected as sources. \n", eq1);
		}
	}


	//----------------------------find completeness magnitude and-------------------------//
	//------------------------------select only events with M>Mc--------------------------//


	if (!cat || (*cat).Mc>=20){

		for (int i=1; i<=eq2; i++){
			eq=seleq2[i-1];
			mag2[i]=mag[eq];
		}
		(*cat).Mc=Mc_maxcurv(mag2+1, eq2)+Mc_offset;

		k=0;
		for (int i=0; i<eq2; i++){
			if (mag2[i+1]>=(*cat).Mc){
				old2new[i+1]=k+1;
				seleq2[k]=seleq2[i];	//doesn't overwrite since k<=i.
				k++;
			}
		}
		eq2=k;

		k=0;
		for (int i=0; i<eq1; i++){
			eq=seleq1[i];
			if (mag[eq]>=(*cat).Mc){
				seleq1[k]=seleq1[i];	//doesn't overwrite since k<=i.
				catindex[k]=old2new[catindex[i]];
				k++;
			}
		}
		eq1=k;

		if(procId == 0) {
			if (flog){
				fprintf(flog, "Calculated completeness magnitude (using maximum curvature): Mc=%.2lf\n", (*cat).Mc);
				fprintf(flog, "%d events selected for LL inversion. \n", eq2);
				fprintf(flog, "%d events selected as sources. \n", eq1);
			}
		}
	}

	//------------------------------fill in catalog:-------------------------//

	if (cat){
		(*cat).pcrst=&crst;
		init_cat1(cat, eq2, gridPMax);

		for (int i=1; i<=eq2; i++){		//todo parallel
		eq=seleq2[i-1];
			SD=f*herr[i];
			SDd=f*verr[i];
			(*cat).t[i]=times[eq];
			(*cat).mag[i]=mag[eq];
			(*cat).lat0[i]=lat[eq];
			(*cat).lon0[i]=lon[eq];
			(*cat).depths0[i]=dep[eq];
			(*cat).err[i]=herr[eq];
			(*cat).verr[i]=verr[eq];
			*((*cat).ngrid + i)=0;
			latlon2localcartesian(lat[eq], lon[eq], crst.lat0, crst.lon0, &y, &x);
			(*cat).x0[i]=x;
			(*cat).y0[i]=y;
	    	if (findgridpoints){
				errP+=find_gridpoints(ygrid, xgrid, dAgrid, depgrid, N, gridPMax, y, x, SD, dep[i], SDd, cut_sd, (*cat).ngrid + i, (*cat).ngridpoints[i], (*cat).weights[i], 1);
				if (errP) break;
			}
		}
		(*cat).tstart=fmax(t0c, (*cat).t[1]);
		(*cat).tend=fmin(t1c, (*cat).t[(*cat).Z]);

		if (!errP) {
		if ((*cat).Mc>=20) (*cat).Mc=Mc_maxcurv((*cat).mag+1, (*cat).Z)+Mc_offset;
		(*cat).b=calculatebvalue((*cat).mag+1, (*cat).Z, (*cat).Mc);

		if(procId == 0) {
			if (verbose_level>1) printf("Estimated GR values for catalog: Mc=%.2lf, b=%.3lf\n", (*cat).Mc, (*cat).b);
				if (flog) {
					fprintf(flog, "Estimated GR values for catalog: Mc=%.2lf, b=%.3lf\n", (*cat).Mc, (*cat).b);
					fflush(flog);
				}
			}
		}
	}

	//------------------------------fill in eqkfm:-------------------------//


	if (eqfm){

		if (Ntot) *Ntot=eq1;
		*eqfm=eqkfm_array(0,eq1-1);

		for (int i=0; i<eq1; i++){	//todo parallel.
			eq=seleq1[i];
			(*eqfm)[i].t=times[eq];
			(*eqfm)[i].lat=lat[eq];
			(*eqfm)[i].lon=lon[eq];
			(*eqfm)[i].depth=dep[eq];
			(*eqfm)[i].mag=mag[eq];
			(*eqfm)[i].index_cat= catindex[i];

	    	latlon2localcartesian((*eqfm)[i].lat, (*eqfm)[i].lon, crst.lat0, crst.lon0, &((*eqfm)[i].y), &((*eqfm)[i].x));
	    	if (findgridpoints){
				if (catindex[i]!=0) errP+=find_gridpoints_d(ygrid, xgrid, depgrid, (*cat).ngridpoints[catindex[i]], (*cat).ngrid[catindex[i]], N, (*eqfm)[i].y, (*eqfm)[i].x, (*eqfm)[i].depth,  (*eqfm)[i].mag, dDCFS,  &((*eqfm)[i].nsel), &((*eqfm)[i].selpoints));
				else errP+=find_gridpoints_d(ygrid, xgrid, depgrid, (int *) 0, 0, N, (*eqfm)[i].y, (*eqfm)[i].x, (*eqfm)[i].depth,  (*eqfm)[i].mag, dDCFS,  &((*eqfm)[i].nsel), &((*eqfm)[i].selpoints));
				if (errP) break;
	    	}

		}
	}

	else {
		if (Ntot) *Ntot=0;
	}

	return (errP!=0);
}


