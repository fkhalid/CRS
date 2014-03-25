/*
 * readcatalog.c
 *
 *  Created on: Dec 22, 2011
 *      Author: camcat
 */

#include "read_RS.h"

int gridPMax=1000;	// max no. points associated with event.

//todo delete this file, is only here for backward compatibility with Parkfield.

int read_RS(char *filename, struct catalog *cat, struct crust crst, double time0, double tstart, double tW1, double tW2, double tend, double Mag, struct eqkfm **eqfm2, double dDCFS, int *nev, int sources){
//sources indicate if events should be used as sources (boolean "aftershocks" in main).
//nev is the total number of events as sources is sources==1, else the total number of events in catalog.

	int err=0, errP;
	int S=countcol(filename), is_in_catalog;
	int inside, *selpts, *sel, eq;
	long Z1,Z0, Zsel;
	double **datas;
	double Lat,Lon,SD, SDlat, SDlon, Depth, SDd, x, y;
	double Latmax, Latmin, Lonmax, Lonmin, Depthmax, Depthmin;
	double toll=1e-10;
	double slip;
	double f=1;	//SD is multiplied by this factor.
	int n_tw=0, Zmax=countline(filename);
	printf("Reading earthquake catalog...\n\n");

	struct eqkfm* eqfm1;
	double *n1, *n2, *sl1, *sl2, *stress1, *stress2;
	double sigma1, sigma2, tau1, tau2;

	if (Zmax<0) return 1;

	//wider areas used to select sources:
	//NB needs to be consistent with readfocmec.c
	Latmax=crst.latmax+0.1;
	Latmin=crst.latmin-0.1;
	Lonmax=crst.lonmax+0.1;
	Lonmin=crst.lonmin-0.1;
	Depthmax=crst.depmax+10;
	Depthmin=fmax(0,crst.depmin-10);

	datas = dmatrix(1, S, 1, Zmax);
	read_matrix(filename, S, 0, datas, &Z1);

	if (Z1>Zmax){
		printf("Wrong number of lines found in catalog %s (Z1=%ld, Zmax=%d). Exiting. \n",filename, Z1, Zmax);
		return (1);
	}

	int eq00=1;
	while (eq00<Z1 && datas[2][eq00]<tstart-toll) eq00++;
	if (eq00==Z1) {
		printf("All events in catalog are occur before starting time. Exiting.\n");
		return (1);
	}

	int eq1=0;	//inside large area.
	int eq2=0;	//inside small area.
	double fr=0.1;
	selpts=ivector(0,Z1);	//contains indices of selected items (obsolete)
	sel=ivector(0,Z1);		//integer sorting selected items (to be used in catalog).
	for (int eq0=eq00; eq0<Z1; eq0++){
		//	if (Lat+SDlat>=cat->latmin && Lat-SDlat<=cat->latmax && Lon+SDlon>=cat->lonmin && Lon-SDlon<=cat->lonmax && Depth+SDd>=cat->depthmin && Depth-SDd<=cat->depthmax){
		//	introduced fr ("fraction") because was getting too many points outside range (Vfrac>1).
		datas[2][eq0]+=time0;
		Lat=datas[3][eq0];
		Lon=datas[4][eq0];
		Depth=datas[5][eq0];
		SD=f*datas[7][eq0];
		SDd=f*datas[8][eq0];
		SDlat=180*SD/(Re*PI);
		SDlon=180*SD/(Re*PI*cos(Lat*pi/180));
		sel[eq0]=0;
		if (datas[6][eq0]>= Mag && datas[2][eq0]<= tend && Lat+fr*SDlat>=Latmin && Lat-fr*SDlat<=Latmax && Lon+fr*SDlon>=Lonmin && Lon-fr*SDlon<=Lonmax && Depth+fr*SDd>=Depthmin && Depth-fr*SDd<=Depthmax){
			selpts[eq1]=eq0;
			eq1+=1;
			if ((datas[2][eq0]<tW1 || datas[2][eq0]> tW2) && Lat+fr*SDlat>=crst.latmin && Lat-fr*SDlat<=crst.latmax && Lon+fr*SDlon>=crst.lonmin && Lon-fr*SDlon<=crst.lonmax && Depth+fr*SDd>=crst.depmin && Depth-fr*SDd<=crst.depmax){
				sel[eq0]=eq2+1;
				eq2+=1;
			}
		}
	}
	Zsel=eq1;
	if (nev) *nev=(sources==1)? Zsel : eq2-n_tw;
	if (sources) eqfm1 = eqkfm_array(0, *nev-1);

	init_cat1(cat, eq2, gridPMax);

	int counter=0;
	#pragma omp parallel for private(counter, errP, x, y, eq, Lat, Lon, Depth, SD, SDd, SDlat, SDlon, is_in_catalog, Z0, inside, slip, n1, n2, sl1, sl2, stress1, stress2, sigma1, sigma2, tau1, tau2)
	for (int eq0=0; eq0<Zsel; eq0++){
		counter+=1;
//		if (omp_get_thread_num()==1) printf("\r %d%% ...     ", (int)(100*(omp_get_num_threads()*counter/Zsel)));
		if (err>0) continue;	//nb: can't break parallel loop.
		errP=0;
		eq=selpts[eq0];
		Lat=datas[3][eq];
		Lon=datas[4][eq];
		Depth=datas[5][eq];
		SD=f*datas[7][eq];
		SDd=f*datas[8][eq];

		SDlat=180*SD/(Re*PI);
		SDlon=180*SD/(Re*PI*cos(Lat*PI/180));

		is_in_catalog=0;
		if (sel[eq]!=0){
			is_in_catalog=1;
			Z0=(datas[2][eq]<tW1)? sel[eq] : sel[eq]-n_tw;		//Z0!=Z1 if some filtering is introduced
			//	if ((datas[2][eq]>0.5 && datas[2][eq]<200) || (datas[2][eq]>-100 && datas[2][eq]<0)) printf("%.0lf\n",datas[1][eq]);
			cat->t[Z0]=datas[2][eq];
			cat->mag[Z0]=datas[6][eq];
			cat->lat0[Z0]=Lat;
			cat->lon0[Z0]=Lon;
			cat->depths0[Z0]=Depth;
			cat->err[Z0]=SD;
			cat->verr[Z0]=SDd;
			*(cat->ngrid + Z0)=0;
			inside =0;
			if (Lat>=crst.latmin && Lat<=crst.latmax && Lon>=crst.lonmin && Lon<=crst.lonmax) inside=1;
			latlon2localcartesian(Lat, Lon, 0.5*(crst.latmax+crst.latmin), 0.5*(crst.lonmax+crst.lonmin), &y, &x);
			cat->x0[Z0]=x;
			cat->y0[Z0]=y;
			errP+=find_gridpoints(crst.y, crst.x, crst.dAgrid, crst.depth, crst.N_allP, gridPMax, y, x, SD, Depth, SDd, 3.0,  cat->ngrid + Z0, cat->ngridpoints[Z0], cat->weights[Z0], inside);
			if (*(cat->ngrid + Z0)==0){
				printf("*** Warning: no grid points selected for event eq=%d! (%lf,%lf,%lf' SD=%lf)\n",eq,Lat,Lon,Depth, SD);
			}
		}

		//TODO make this function more general and use it even for the case aftershocks!=1?

		if (eqfm2){
			if (sources==1){
				eqfm1[eq0].t=datas[2][eq];
				eqfm1[eq0].lat=datas[3][eq];
				eqfm1[eq0].lon=datas[4][eq];
				eqfm1[eq0].depth=datas[5][eq];
				eqfm1[eq0].mag=datas[6][eq];
				eqfm1[eq0].noise=0;
				eqfm1[eq0].index_cat= (is_in_catalog==1)? Z0 : 0;
				latlon2localcartesian(eqfm1[eq0].lat, eqfm1[eq0].lat, 0.5*(crst.latmax+crst.latmin), 0.5*(crst.lonmax+crst.lonmin), &(eqfm1[eq0].y), &(eqfm1[eq0].x));
				if (is_in_catalog==1) errP+=find_gridpoints_d(crst.y, crst.x, crst.depth, cat->ngridpoints[Z0], cat->ngrid[Z0], crst.N_allP, eqfm1[eq0].y, eqfm1[eq0].x, eqfm1[eq0].depth,  eqfm1[eq0].mag, dDCFS,  &(eqfm1[eq0].nsel), &(eqfm1[eq0].selpoints));
				else errP+=find_gridpoints_d(crst.y, crst.x, crst.depth, (int *) 0, 0, crst.N_allP, eqfm1[eq0].y, eqfm1[eq0].x, eqfm1[eq0].depth,  eqfm1[eq0].mag, dDCFS,  &(eqfm1[eq0].nsel), &(eqfm1[eq0].selpoints));
				switch (S){

				case 14:
					eqfm1[eq0].is_slipmodel=1;
					eqfm1[eq0].np_st=1;
					eqfm1[eq0].np_di=1;
					eqfm1[eq0].slip_str=dvector(1,1);
					eqfm1[eq0].slip_dip=dvector(1,1);
					eqfm1[eq0].pos_s=dvector(1,1);	//location of patches within fault; [0], [0] for single patch events.
					eqfm1[eq0].pos_d=dvector(1,1);
					eqfm1[eq0].pos_s[1]=0;	//location of patches within fault; [0], [0] for single patch events.
					eqfm1[eq0].pos_d[1]=0;
					eqfm1[eq0].str1=datas[9][eq];
					eqfm1[eq0].dip1=datas[10][eq];
					eqfm1[eq0].rake1=fmod(datas[11][eq]+360.0,360.0);
					eqfm1[eq0].str2=datas[12][eq];
					eqfm1[eq0].dip2=datas[13][eq];
					eqfm1[eq0].rake2=fmod(datas[14][eq]+360.0,360.0);	//so don't have -ve values.

					//TODO in theory, should find L and W for both mech. (rake differs).
					WellsCoppersmith(eqfm1[eq0].mag, eqfm1[eq0].rake1, &(eqfm1[eq0].L), &(eqfm1[eq0].W), &slip);
					eqfm1[eq0].tot_slip=pow(10,(1.5*(eqfm1[eq0].mag+6)))*(1.0/(crst.mu*pow(10,12)*eqfm1[eq0].W*eqfm1[eq0].L));
					eqfm1[eq0].whichfm=0;	//this means to preference, hence don't calculate eqfm1[eq].slip_str[1], eqfm1[eq].slip_dip[1]. //TODO need to make this consistent in future functions.

					//find which plane is best oriented to regional stress field (this part is a bottleneck when repeated over thousands of f.m.):
					//n1[1]=-sin(eqfm1[eq].str1*DEG2RAD)*sin(eqfm1[eq].dip1*DEG2RAD);
					//n1[2]=cos(eqfm1[eq].str1*DEG2RAD)*sin(eqfm1[eq].dip1*DEG2RAD);
					//n1[3]=-cos(eqfm1[eq].dip1*DEG2RAD);
					//n2[1]=-sin(eqfm1[eq].str2*DEG2RAD)*sin(eqfm1[eq].dip2*DEG2RAD);
					//n2[2]=cos(eqfm1[eq].str2*DEG2RAD)*sin(eqfm1[eq].dip2*DEG2RAD);
					//n2[3]=-cos(eqfm1[eq].dip2*DEG2RAD);
					//
					//sl1[1]=cos(eqfm1[eq].str1*DEG2RAD)*cos(eqfm1[eq].rake1*DEG2RAD)+sin(eqfm1[eq].str1*DEG2RAD)*sin(eqfm1[eq].rake1*DEG2RAD)*cos(eqfm1[eq].dip1*DEG2RAD);
					//sl1[2]=sin(eqfm1[eq].str1*DEG2RAD)*cos(eqfm1[eq].rake1*DEG2RAD)-cos(eqfm1[eq].str1*DEG2RAD)*sin(eqfm1[eq].rake1*DEG2RAD)*cos(eqfm1[eq].dip1*DEG2RAD);
					//sl1[3]=-sin(eqfm1[eq].rake1*DEG2RAD)*sin(eqfm1[eq].dip1*DEG2RAD);
					//sl2[1]=cos(eqfm1[eq].str2*DEG2RAD)*cos(eqfm1[eq].rake2*DEG2RAD)+sin(eqfm1[eq].str2*DEG2RAD)*sin(eqfm1[eq].rake2*DEG2RAD)*cos(eqfm1[eq].dip2*DEG2RAD);
					//sl2[2]=sin(eqfm1[eq].str2*DEG2RAD)*cos(eqfm1[eq].rake2*DEG2RAD)-cos(eqfm1[eq].str2*DEG2RAD)*sin(eqfm1[eq].rake2*DEG2RAD)*cos(eqfm1[eq].dip2*DEG2RAD);
					//sl2[3]=-sin(eqfm1[eq].rake2*DEG2RAD)*sin(eqfm1[eq].dip2*DEG2RAD);
					//
					//mtimesv(crst.S,n1,stress1,3,3);
					//mtimesv(crst.S,n2,stress2,3,3);
					//vdotv(stress1,n1,&sigma1,3);
					//vdotv(stress2,n2,&sigma2,3);
					//vdotv(stress1,sl1,&tau1,3);
					//vdotv(stress2,sl2,&tau2,3);
					//
					//if (tau1+crst.fric*sigma1>tau2+crst.fric*sigma2){
					//eqfm1[eq0].slip_str[1]=slip*cos(DEG2RAD*eqfm1[eq0].rake1);
					//eqfm1[eq0].slip_dip[1]=slip*sin(DEG2RAD*eqfm1[eq0].rake1);
					//eqfm1[eq0].whichfm=1;
					//}
					//else{
					//	eqfm1[eq0].slip_str[1]=slip*cos(DEG2RAD*eqfm1[eq0].rake2);
					//	eqfm1[eq0].slip_dip[1]=slip*sin(DEG2RAD*eqfm1[eq0].rake2);
					//	eqfm1[eq0].whichfm=2;
					//}
					break;

				case 11:
					eqfm1[eq0].is_slipmodel=1;
					eqfm1[eq0].np_st=1;
					eqfm1[eq0].np_di=1;
					eqfm1[eq0].slip_str=dvector(1,1);
					eqfm1[eq0].slip_dip=dvector(1,1);
					eqfm1[eq0].pos_s=dvector(1,1);	//location of patches within fault; [0], [0] for single patch events.
					eqfm1[eq0].pos_d=dvector(1,1);
					eqfm1[eq0].pos_s[1]=0;	//location of patches within fault; [0], [0] for single patch events.
					eqfm1[eq0].pos_d[1]=0;
					eqfm1[eq0].str1=datas[9][eq];
					eqfm1[eq0].dip1=datas[10][eq];
					eqfm1[eq0].rake1=fmod(datas[11][eq]+360.0,360.0);
					WellsCoppersmith(eqfm1[eq0].mag, eqfm1[eq0].rake1, &(eqfm1[eq0].L), &(eqfm1[eq0].W), &slip);
					eqfm1[eq0].tot_slip=pow(10,(1.5*(eqfm1[eq0].mag+6)))*(1.0/(crst.mu*pow(10,12)*eqfm1[eq0].W*eqfm1[eq0].L));
					eqfm1[eq0].whichfm=1;
					eqfm1[eq0].slip_str[1]=slip*cos(DEG2RAD*eqfm1[eq0].rake1);
					eqfm1[eq0].slip_dip[1]=slip*sin(DEG2RAD*eqfm1[eq0].rake1);
					break;

				case 8:
					eqfm1[eq0].is_slipmodel=0;
					eqfm1[eq0].np_st= eqfm1[eq0].np_di=0;
					break;

				default:
					printf("*Wrong number of columns (%d) in input file: %s*\n", S, filename);
				}
			}

			else {
				//this is only used to calculate rates in forecast function.
				if (is_in_catalog==1){
					eqfm1[Z0-1].is_slipmodel=0;
					eqfm1[Z0-1].t=datas[2][eq];
					eqfm1[Z0-1].lat=datas[3][eq];
					eqfm1[Z0-1].lon=datas[4][eq];
					eqfm1[Z0-1].depth=datas[5][eq];
					eqfm1[Z0-1].mag=datas[6][eq];
					eqfm1[Z0-1].noise=0;
					eqfm1[Z0-1].index_cat= Z0;
					eqfm1[Z0-1].selpoints=cat->ngridpoints[Z0];
					eqfm1[Z0-1].nsel=cat->ngrid[Z0];
				}
			}
		}

	if (errP>0) {
			#pragma omp atomic
			err+=1;
		}
	}

	if (err>0) {printf("*** %d errors in parallel section of readcatalogandfocmec.c. Exiting. ***\n",err); return(1);}
	cat->tstart=cat->t[1];
	cat->tend=tend;

	if (eqfm2)*eqfm2=eqfm1;
	free_dmatrix(datas, 1, S, 1, Zmax);
	free_ivector(selpts, 0, Z1);
	free_ivector(sel, 0, Z1);
	printf("done.\n %ld events found, %ld used in catalog,  %ld used as sources.\n\n",Z1, cat->Z, sources*Zsel);

	return(0);

}
