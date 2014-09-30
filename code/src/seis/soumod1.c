// soumod.c reads slip distribution with a certain discretisation step from given file sourceold, 
// resamples it to the given discretisation step disc and adjust the values that way 
// that a k²-distribution of slip is obtained. Hence, the low frequency content of the original
// slip distribution is kept and the high frequency content is varied randomly.
//
// if no slip distribution is given, the programm creates a random slip distribution with the parameters given
//
// The focal mechanisms of the different fault patches are varied. Variation of strike and dip are chosen, so that 
// the fault plane assumes the shape of a fractal surface with aspect ratio 0.1 (in respect to the shorter fault edge)
// Rake is varied by +-10 degree
//
//////////////////////////////////////////////// The rupture velocity is varied.
//
// Katrin Kieling, august 2010 
// modified by Camilla Cattania, March 2013
//
// last modified: 03/08/2010 (Katrin), 13/03/2013 (Camilla).
// 

#include "soumod1.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

double tot_slip (struct eqkfm eqfm1){

	double slip=0.0;
	int NP=eqfm1.np_di*eqfm1.np_st;

	for (int p=1; p<=NP; p++) {
		slip+=pow(eqfm1.slip_str[p]*eqfm1.slip_str[p]+eqfm1.slip_dip[p]*eqfm1.slip_dip[p],0.5);
	}

	return slip/NP;
}

int scale_to_mag(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double * slips, double *rakes){

	double slip, M0old, M0;	//slip, seismic moments.
	int N1= eqkfm1.np_di*eqkfm1.np_st;
	int N2= (*eqkfm2).np_di*(*eqkfm2).np_st;

	M0old=0.0;
	for (int p=1; p<=N1; p++) {
		slip=pow(eqkfm1.slip_str[p]*eqkfm1.slip_str[p]+eqkfm1.slip_dip[p]*eqkfm1.slip_dip[p],0.5);
		M0old+=slip/N1;
	}

	M0=0.0;
	for (int p=1; p<=N2; p++) M0+=slips[p]/N2;

	for (int p=1; p<=N2; p++) {
		slip=slips[p]*M0old/M0;
		(*eqkfm2).slip_str[p]=slip*cos((rakes[p])*pi/180);
		(*eqkfm2).slip_dip[p]=-slip*sin((rakes[p])*pi/180);
		if (isnan((*eqkfm2).slip_dip[p]) || isnan((*eqkfm2).slip_str[p])) return 1;
	}

	return 0;
}

int adjust_faults(struct eqkfm *eqkfm0,  int NF, int vert){
	//vert=1: vertical; vert=0: horizontal.

	double *dipV, *P;
	double newlat, newlon, dx, dy, dip, str, diff;

	dipV=dvector(1,3);
	P=dvector(1,3);

	for (int n=1; n<=NF-1; n++){

		dip=eqkfm0[n-1].dip1*DEG2RAD;
		str=eqkfm0[n-1].str1*DEG2RAD;

		dipV[1]=cos(dip)*sin(str+PI/2.0);
		dipV[2]=cos(dip)*cos(str+PI/2);
		dipV[3]=sin(dip);

		for (int i=1; i<=3; i++) P[i]=dipV[i]*eqkfm0[n-1].W;
		P[3]+=eqkfm0[n-1].depth;
		localcartesian2latlon(P[2],P[1],eqkfm0[n-1].lat,eqkfm0[n-1].lon, &newlat, &newlon);
		latlon2localcartesian(eqkfm0[n].lat,eqkfm0[n].lon, newlat, newlon, &dy, &dx);
		diff=sqrt(pow(dx,2)+pow(dy,2)+pow(eqkfm0[n].depth-P[3],2));

		if (extra_verbose) {
			print_screen("Fault adjustment (fault no.%d): [%lf,%lf,%lf]  -> [%lf,%lf,%lf]\n", n, eqkfm0[n].lat,eqkfm0[n].lon, eqkfm0[n].depth, newlat, newlon, P[3]);
			print_logfile("Fault adjustment (fault no.%d): [%lf,%lf,%lf]  -> [%lf,%lf,%lf]\n", n, eqkfm0[n].lat,eqkfm0[n].lon, eqkfm0[n].depth, newlat, newlon, P[3]);
		}
		if (diff>1.0) {
			if (diff>10.0){
			print_screen("*** Error - shifting subfault no.%d by more than 5km in adjust_faults. *** \n",n);
			print_logfile("*** Error - shifting subfault no.%d by more than 5km in adjust_faults. *** \n",n);
			return (1);
			}
			print_screen("*** Warning - shifting subfault no.%d by more than 1km in adjust_faults. *** \n",n);
			print_logfile("*** Warning - shifting subfault no.%d by more than 1km in adjust_faults. *** \n",n);
		}

		eqkfm0[n].lat=newlat;
		eqkfm0[n].lon=newlon;
		eqkfm0[n].depth=P[3];
	}

	free_dvector(dipV,1,3);
	free_dvector(P,1,3);

	return(0);
}

int which_taper(struct eqkfm *eqkfm0,  int NF, int tap_bot, int tap_top, double d_close_enough){
//setup taper vectors in eqkfm0 (taper [1...4]=top,bottom,right,left) depending on how faults are oriented relative to each others.

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double *P;
	double cosH, H, Hdip;
	int Ver_align=0, Hor_align=0;

	if (NF==1){
		if (eqkfm0[0].taper==NULL) eqkfm0[0].taper=ivector(1,4);
		for (int i=1; i<=4; i++) eqkfm0[0].taper[i]=1;
		return(0);
	}

	else {
		P=dvector(1,3);

		for (int n=2; n<=NF; n++) {
			latlon2localcartesian(eqkfm0[n-1].lat,eqkfm0[n-1].lon, eqkfm0[n-2].lat, eqkfm0[n-2].lon, P+2, P+1);
			P[3]=eqkfm0[n-1].depth-eqkfm0[n-2].depth;
			cosH=P[3]/sqrt(pow(P[1],2)+pow(P[2],2));
			H=atan(cosH)*RAD2DEG;
			Hdip=eqkfm0[n-2].dip1;
			if (H<1.0) Hor_align+=1;
			else {
				if (fabs(H-Hdip)<1.0) Ver_align+=1;
			}
		}

		for (int n=1; n<NF-1; n++){
			if (eqkfm0[n].taper==NULL) eqkfm0[n].taper=ivector(1,4);
		}

		if (Hor_align==NF-1 && Ver_align==0){

			print_screen("Mainshock faults oriented horizontally.\n");
			print_logfile("Mainshock faults oriented horizontally.\n");

			for (int n=0; n<NF; n++){
				eqkfm0[n].taper[4]=1;
				eqkfm0[n].taper[3]=1;
				eqkfm0[n].taper[1]=tap_top;
				eqkfm0[n].taper[2]=tap_bot;
			}

			double d, d0=1e10;
			double east, north;
			int i0, j0;

			for (int n=2; n<=NF; n++) {
				d0=1e10;
				latlon2localcartesian(eqkfm0[n-1].lat,eqkfm0[n-1].lon, eqkfm0[n-2].lat, eqkfm0[n-2].lon, P+2, P+1);
				for (int i=-1; i<=1; i+=2){
					for (int j=-1; j<=1; j+=2){
						north= 0.5*(i*eqkfm0[n-1].L*cos(DEG2RAD*eqkfm0[n-1].str1)+j*eqkfm0[n-2].L*cos(DEG2RAD*eqkfm0[n-2].str1));
						east= 0.5*(i*eqkfm0[n-1].L*sin(DEG2RAD*eqkfm0[n-1].str1)+j*eqkfm0[n-2].L*sin(DEG2RAD*eqkfm0[n-2].str1));
						d=pow((north-P[2])*(north-P[2])+(east-P[1])*(east-P[1]),0.5);
						if (d<d0){
							d0=d;
							i0=i;
							j0=j;
						}
					}
				}

				if (d0<=d_close_enough){	//only avoid tapering if the faults are close enough.

					if (j0==1) eqkfm0[n-2].taper[3]=0;		//right;
					else eqkfm0[n-2].taper[4]=0;			//left;

					if (i0==1) eqkfm0[n-1].taper[4]=0;		//left;
					else eqkfm0[n-1].taper[3]=0;			//right;
				}
			}

	//		if (sign>0){
	//			eqkfm0[0].taper[4]=1;
	//			eqkfm0[0].taper[3]=0;
	//			eqkfm0[NF-1].taper[3]=1;
	//			eqkfm0[NF-1].taper[4]=0;
	//		}
	//		else {
	//			eqkfm0[0].taper[4]=0;
	//			eqkfm0[0].taper[3]=1;
	//			eqkfm0[NF-1].taper[3]=0;
	//			eqkfm0[NF-1].taper[4]=1;
	//		}
	//
	//		eqkfm0[0].taper[1]=tap_top;
	//		eqkfm0[0].taper[2]=tap_bot;
	//		eqkfm0[NF-1].taper[1]=tap_top;
	//		eqkfm0[NF-1].taper[2]=tap_bot;
		}

		else {
			if (Ver_align==NF-1 && Hor_align==0){

				print_screen("Mainshock faults oriented vertically.\n");
				print_logfile("Mainshock faults oriented vertically.\n");
				adjust_faults(eqkfm0,  NF, 1);

				for (int n=1; n<NF-1; n++){
					eqkfm0[n].taper[1]=0;
					eqkfm0[n].taper[2]=0;
					eqkfm0[n].taper[3]=1;
					eqkfm0[n].taper[4]=1;
				}
				eqkfm0[0].taper[1]=tap_top;
				eqkfm0[NF-1].taper[2]=tap_bot;
			}

			else {
				print_screen("** Warning: could not determine relative orientation of faults - model is not tapered! **\n");
				print_logfile("** Warning: could not determine relative orientation of faults - model is not tapered! **\n");
				for (int n=0; n<NF-1; n++){
					for (int i=1; i<=4; i++) eqkfm0[n].taper[i]=0;
				}
			}

		}

		free_dvector(P,1,3);
		return(0);
	}

}

int suomod1_resample(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double disc, double velmean){
  char     outname2[]="source_res.out";
  double   reflat, reflon, refdepth, strike, rake, dip1;
  double   odiscx, odiscy, maxox=-1e10, maxoy=-1e10, minox=1e10, minoy=1e10, start_x, start_y;
  double   *ox, *oy, *REslipo, *rake_v, *strike_v, *dip;
  double   *nx, *ny, *REold, *IMold, *dipold, *rakeold, *strikeold;
  double   oxdim=0.0 /*, oxdim_tot*/;
  double   oydim=0.0;
  double 	*slip;
  double 	ndiscx, ndiscy;
  int	    print_Fourier, printout;
  int      i, j, k, l;
  int      ns, nns, nsx, nsy;

  //can be set manually to print out extra files.
  printout= 0;
  print_Fourier= 0;

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

    ns=eqkfm1.np_st*eqkfm1.np_di;
    REslipo= dvector(1,ns);
    rake_v= dvector(1,ns);
    strike_v= dvector(1,ns);
    dip= dvector(1,ns);
	reflat=eqkfm1.lat;
	reflon=eqkfm1.lon;
	refdepth=eqkfm1.depth;
	switch (eqkfm1.whichfm){
		case 1:
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip1=eqkfm1.dip1;
			break;
		case 2:
			rake=eqkfm1.rake2;
			strike=eqkfm1.str2;
			dip1=eqkfm1.dip2;
			break;
		case 0:
			print_screen("Warning: ambiguous focal plane for eqkfm1 (suomod1_resample) -> using first plane!\n");
			print_logfile("Warning: ambiguous focal plane for eqkfm1 (suomod1_resample) -> using first plane!\n");
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip1=eqkfm1.dip1;
			break;
		default:
			print_screen("Error: eqkfm1.whichfm has illegal value in suomod1_resample(%d)!\n",eqkfm1.whichfm);
			print_logfile("Error: eqkfm1.whichfm has illegal value in suomod1_resample(%d)!\n",eqkfm1.whichfm);
			return(1);
	}
	oxdim=eqkfm1.L;
	oydim=eqkfm1.W;
	nsx=eqkfm1.np_st;
	nsy=eqkfm1.np_di;
	ox=eqkfm1.pos_s;
	oy=eqkfm1.pos_d;
	odiscx=oxdim/nsx;
	odiscy=oydim/nsy;

	for (int np=1; np<=ns; np++){
		REslipo[np]=pow(eqkfm1.slip_str[np]*eqkfm1.slip_str[np]+eqkfm1.slip_dip[np]*eqkfm1.slip_dip[np],0.5);
		rake_v[np]=(-180/pi)*atan2(eqkfm1.slip_dip[np],eqkfm1.slip_str[np]);	//check this! (sign)
		if (isnan(rake_v[np])==1) rake_v[np]=rake;	//for example if slip=0 in both directions.
		strike_v[np]=strike;
		dip[np]=dip1;
	}

	for (i=1;i<=ns;i++){
	   if (ox[i]>maxox) maxox=ox[i];
	   if (oy[i]>maxoy) maxoy=oy[i];
	   if (ox[i]<minox) minox=ox[i];
	   if (oy[i]<minoy) minoy=oy[i];
	}

	//oxdim_tot=(maxox-minox)+odiscx;	//could be larger than oxdim if fault is skewed. Is this needed?

	 maxoy=maxoy+odiscy;
	 //resampling of slipmap
	 //todo should use oxdim or oxdim_tot?

	  nns=(int)(ceil(oxdim/disc)*ceil(oydim/disc));
	  nsx=(int)ceil(oxdim/disc);
	  nsy=(int)ceil(oydim/disc);
	  ndiscx=oxdim/nsx;
	  ndiscy=oydim/nsy;
	  start_x=minox-0.5*odiscx+0.5*ndiscx;
	  start_y=minoy-0.5*odiscy+0.5*ndiscy;
	  REold    = dvector(1,nns);
	  IMold    = dvector(1,nns);
	  dipold   = dvector(1,nns);
	  rakeold   = dvector(1,nns);
	  strikeold   = dvector(1,nns);
	  nx       = dvector(1,nns);
	  ny       = dvector(1,nns);
      slip=dvector(1,nns);

      	  i=1;
	  for (k=1; k<=nsy; k++)
	        { for (j=1;j<=nsx;j++)
	                { nx[i]=start_x+(j-1)*ndiscx; //relative position along strike (note: new fault not skewed).
	                  ny[i]=start_y+(k-1)*ndiscy; //relative position along dip
	                  i++;
	                }
	        }
      if (ndiscx<odiscx && ndiscy<odiscy) {// refine slip distribution
          for (int i=1; i<=nns; i++){
        	  l=1;
        	  while (l<ns && (fabs(ox[l]-nx[i])>odiscx || fabs(oy[l]-ny[i])>odiscy)) l++;
        	  if (l==ns && (fabs(ox[l]-nx[i])>odiscx || fabs(oy[l]-ny[i])>odiscy)) {	//new point is outside old fault (possible for skewed fault).
            	  REold[i]=0.0;
    			  IMold[i]=0;
    			  dipold[i]=dip1;
    			  rakeold[i]=rake;
    			  strikeold[i]=strike;
        	  }
        	  else {
				  REold[i]=REslipo[l];	//new point is inside new fault
				  IMold[i]=0;
				  dipold[i]=dip[l];
				  rakeold[i]=rake_v[l];
				  strikeold[i]=strike_v[l];
        	  }
          }
      }

	  else if (ndiscx>odiscx && ndiscy>odiscy){
		  if (extra_verbose) {
			  print_screen("New resolution is larger than old one: will not resample (suomod1_resample).\n");
			  print_logfile("New resolution is larger than old one: will not resample (suomod1_resample).\n");
		  }
		  copy_eqkfm_all(eqkfm1, eqkfm2);
		  return (0);
	  }

	  else  {
		  if (extra_verbose) {
			  print_screen("Model has right discretization - will not be resampled. \n");
		  }
		  copy_eqkfm_all(eqkfm1, eqkfm2);
		  return(0);
	  }

  //--------------Fill in eqkfm2.-------------------//

    copy_eqkfm_noslipmodel(eqkfm1, eqkfm2);

 	(*eqkfm2).L=oxdim;
 	(*eqkfm2).W=oydim;
	(*eqkfm2).np_st=nsx;
	(*eqkfm2).np_di=nsy;
	(*eqkfm2).pos_s=dvector(1,nns);
	(*eqkfm2).pos_d=dvector(1,nns);
	(*eqkfm2).slip_str=dvector(1,nns);
	(*eqkfm2).slip_dip=dvector(1,nns);

	for (int p=1; p<=nns; p++) {
		(*eqkfm2).pos_s[p]=nx[p];
		(*eqkfm2).pos_d[p]=ny[p];
	}

	scale_to_mag(eqkfm1, eqkfm2, REold, rakeold);		//final slip.
	(*eqkfm2).tot_slip=tot_slip(*eqkfm2);


  //--------------Print things out-------------------//

  if (printout==1) {
	  if(procId == 0) {
		  print_slipmodel(outname2,eqkfm2,1);
	  }
  }

  //----------------Free memory----------------------//

  free_dvector(REold,1,nns);
  free_dvector(REslipo,1,ns);
  free_dvector(rake_v,1,ns);
  free_dvector(strike_v,1,ns);
  free_dvector(dip,1,ns);
  free_dvector(IMold,1,nns);
  free_dvector(strikeold,1,nns);
  free_dvector(dipold,1,nns);
  free_dvector(rakeold,1,nns);
  free_dvector(slip,1,nns);
  free_dvector(nx,1,nns);
  free_dvector(ny,1,nns);
//  if (ndiscx>odiscx && ndiscy>odiscy){
//	  free_dvector(lat,1,nns);
//	  free_dvector(lon,1,nns);
//	  free_dvector(z,1,nns);
//  }

  return(0);
}

int suomod1_taper(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
  char     outname2[200];
  double   strike, rake, dip;
  double   odiscx, odiscy, maxox=-1e10, maxoy=-1e10, minox=1e10, minoy=1e10, *start_x, start_y;
  double   *ox, *oy, *REslipo, *slip, *rakes;
  double   *nx, *ny;
  double 	*taper;
  double   alphax, alphay;
  double   oxdim=0.0;
  double   oydim=0.0;
  int	   printout;
  int      ns, nns, nsx, nsy;
  int 	   top, bottom, right, left;	//where to taper.

    // can be set to 1 to print extra files.
	printout=0;

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	top=eqkfm1.taper[1];
	bottom=eqkfm1.taper[2];
	right=eqkfm1.taper[3];
	left=eqkfm1.taper[4];

    ns=nns=eqkfm1.np_st*eqkfm1.np_di;
    if (eqkfm1.np_st==1) left=right=0;
    if (eqkfm1.np_di==1) top=bottom=0;
    if (ns==1) {
    	copy_eqkfm_all(eqkfm1, eqkfm2);	//bug fix: using copy_eqkfm_slipmodel did not copy nsel, so no points were selected.
    	print_screen("** Warning: model has a single patch, will not be tapered.**\n");
    	print_logfile("** Warning: model has a single patch, will not be tapered.**\n");
    	return 0;
    }

	REslipo=dvector(1,ns);
	rakes=dvector(1,ns);
	taper=dvector(1,ns);
	start_x = dvector(1,ns);
	nx = dvector(1,nns);
	ny = dvector(1,nns);
	slip = dvector(1,nns);

    if (printout) sprintf(outname2, "source_tap%d%d%d%d.out",top, bottom, right, left);

	switch (eqkfm1.whichfm){
		case 1:
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip=eqkfm1.dip1;
			break;
		case 2:
			rake=eqkfm1.rake2;
			strike=eqkfm1.str2;
			dip=eqkfm1.dip2;
			break;
		case 0:
			print_screen("Warning: ambiguous focal plane for eqkfm1 (suomod1_taper) -> using first plane!\n");
			print_logfile("Warning: ambiguous focal plane for eqkfm1 (suomod1_taper) -> using first plane!\n");
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip=eqkfm1.dip1;
			break;
		default:
			print_screen("Error: eqkfm1.whichfm has illegal value in suomod1_taper(%d)!\n",eqkfm1.whichfm);
			print_logfile("Error: eqkfm1.whichfm has illegal value in suomod1_taper(%d)!\n",eqkfm1.whichfm);
			return(1);
			break;
	}
	oxdim=eqkfm1.L;
	oydim=eqkfm1.W;
	nsx=eqkfm1.np_st;
	nsy=eqkfm1.np_di;
	ox=eqkfm1.pos_s;
	oy=eqkfm1.pos_d;
	odiscx=oxdim/nsx;
	odiscy=oydim/nsy;

	for (int np=1; np<=ns; np++) {
		REslipo[np]=pow(eqkfm1.slip_str[np]*eqkfm1.slip_str[np]+eqkfm1.slip_dip[np]*eqkfm1.slip_dip[np],0.5);
		rakes[np]=(-180/pi)*atan2(eqkfm1.slip_dip[np],eqkfm1.slip_str[np]);
	}

	for (int i=1;i<=ns;i++){
	   if (ox[i]>maxox) maxox=ox[i];
	   if (oy[i]>maxoy) {
		   maxoy=oy[i];
		   start_x[i]=ox[i];	//new row of patches -> new starting point (used for tapering).
	   }
	   else start_x[i]=start_x[i-1];
	   if (ox[i]<minox) minox=ox[i];
	   if (oy[i]<minoy) minoy=oy[i];
	}
	start_y=minoy;	//assume this is the same for all columns (unlike for rows).
	for (int i=1; i<=ns; i++){
	  nx[i]=ox[i]-start_x[i];
	  ny[i]=oy[i]-start_y;
	}

	//-------Tapering-----------//
	/* create taper (similar to tukey window but with onlay 1/4 of the cos period, alpha = fraction of length used for decay, alpha=1 = rectangular window, alpha=0 = similar Hann window)
	 */

	alphax=0.51;
	alphay=0.51;

	for (int i=1;i<=nns;i++){
		taper[i] = ((right==1 && (nx[i]-(nsx-1)/2*odiscx)>alphax*nsx/2*odiscx) || (left==1 && -(nx[i]-(nsx-1)/2*odiscx)>alphax*nsx/2*odiscx))?
			  cos(pi/2*(fabs(fabs(nx[i]-(nsx-1)/2*odiscx)-alphax*nsx/2*odiscx))/((1-alphax)*nsx/2*odiscx)) : 1;
		if ((bottom==1 && (ny[i]-(nsy-1)/2*odiscy)>alphay*nsy/2*odiscy) || (top==1 && -(ny[i]-(nsy-1)/2*odiscy)>alphay*nsy/2*odiscy)){
			  taper[i] *=(cos(pi/2*(fabs(fabs(ny[i]-(nsy-1)/2*odiscy)-alphay*nsy/2*odiscy))/((1-alphay)*nsy/2*odiscy)));
		}
	}

  //Fill in eqkfm2.
  if (eqkfm2!=&eqkfm1){

	copy_eqkfm_noslipmodel(eqkfm1, eqkfm2);

	(*eqkfm2).L=oxdim;
	(*eqkfm2).W=oydim;
	(*eqkfm2).np_st=nsx;
	(*eqkfm2).np_di=nsy;
	(*eqkfm2).pos_s=eqkfm1.pos_s;
	(*eqkfm2).pos_d=eqkfm1.pos_d;
	if ((*eqkfm2).slip_str==NULL) (*eqkfm2).slip_str=dvector(1,nns);
	if ((*eqkfm2).slip_dip==NULL) (*eqkfm2).slip_dip=dvector(1,nns);

  }

	for (int p=1; p<=nns; p++) slip[p]=REslipo[p]*taper[p];		//final slip.

	scale_to_mag(eqkfm1, eqkfm2, slip, rakes);		//final slip.
	(*eqkfm2).tot_slip=tot_slip(*eqkfm2);

  //--------------Print things out-------------------//


	if (printout==1) {
		if(procId == 0) {
			print_slipmodel(outname2,eqkfm2,1);
		}
	}


  //--------------Free memory-------------------//

  free_dvector(slip,1,ns);
  free_dvector(REslipo,1,ns);
  free_dvector(start_x,1,ns);
  free_dvector(nx,1,ns);
  free_dvector(ny,1,ns);
  free_dvector(taper,1,ns);
  free_dvector(rakes,1,ns);
  return 0;

}

int suomod1_hf(struct eqkfm eqkfm0, struct eqkfm *eqkfm2, double H, long *seed, int noise_only){
/* if noise_only==1, returns noise (not old model + noise). Keeps info about magnitude into new model, which can be used to rescale the noise to future models.
 * if (*eqkfm2).slip_str(dip) are not null, function assumes that they have been previously initialized to correct size, and does not allocate memory. *
 */

	FILE     *fout;
	char 	fname[120];
	double velmean=0.0;
	struct   eqkfm *eqkfm1;
	char     outname2[120];
	float    Me= (float) (eqkfm0.mag);
	double 	toll=1e-10;
	double   strike, rake, dip1;
	double   minslip0, maxslip0;
	double   odiscx, odiscy, maxox=-1e10, maxoy=-1e10, minox=1e10, minoy=1e10, start_x, start_y;
	double   kny_x, kny_y, dk_x, dk_y, kc, kcx, kcy;
	double   k_i, A;
	double   *ox, *oy;
	double   *nx, *ny, *REslip, *IMslip, *k_x, *k_y, *REold, *IMold, *FToldRE, *FToldIM, *dipold, *rakeold, *strikeold;
	double   *REsurf, *IMsurf, *FTsurfRE, *FTsurfIM, *dstrike, *ddip, *drake, *FTrand2RE, *FTrand2IM;
	double   *rakes;
	double   height, maxsurf, minsurf;
	//double   xhypo, yhypo;
	double   *FTrandRE, *FTrandIM, *specRE, *specIM /* ,*vel*/;
	double   minslip, diffspec;
	double   oxdim=0.0;
	double   oydim=0.0;
	double   disc;
	double   aspect=0.1;
	double 	 /*M0, M0old,*/ *slip;
	int	   	 old=1, print_Fourier, printout;
	int      i, j, k, l;
	int      ns, nns, nsx, nsy;
	int	   	 rough=0;
	long     seed2 = *seed;

	//can be set to 1 to print extra output files.
	printout= 0;
	print_Fourier= 0;

  	switch (eqkfm0.whichfm){
		case 1:
			rake=eqkfm0.rake1;
			strike=eqkfm0.str1;
			dip1=eqkfm0.dip1;
			break;
		case 2:
			rake=eqkfm0.rake2;
			strike=eqkfm0.str2;
			dip1=eqkfm0.dip2;
			break;
		case 0:
			print_screen("Warning: ambiguous focal plane for eqkfm1 (suomod1_hf) -> using first plane!\n");
			print_logfile("Warning: ambiguous focal plane for eqkfm1 (suomod1_hf) -> using first plane!\n");
			rake=eqkfm0.rake1;
			strike=eqkfm0.str1;
			dip1=eqkfm0.dip1;
			break;
		default:
			print_screen("Error: eqkfm1.whichfm has illegal value in suomod1_hf(%d)!\n",eqkfm0.whichfm);
			print_logfile("Error: eqkfm1.whichfm has illegal value in suomod1_hf(%d)!\n",eqkfm0.whichfm);
			return(1);
  	}

  	odiscx=eqkfm0.L/eqkfm0.np_st;
  	odiscy=eqkfm0.W/eqkfm0.np_di;
  	if (fabs(odiscx-odiscy)>toll) {
  		eqkfm1=eqkfm_array(0,0);
  		suomod1_resample(eqkfm0, eqkfm1, fmin(odiscx, odiscy), velmean);
  	}
	else eqkfm1=&eqkfm0;

	oxdim=(*eqkfm1).L;
	oydim=(*eqkfm1).W;
	nsx=(*eqkfm1).np_st;
	nsy=(*eqkfm1).np_di;
	ox=(*eqkfm1).pos_s;
	oy=(*eqkfm1).pos_d;
	nns=ns=(*eqkfm1).np_di*(*eqkfm1).np_st;

	disc=oxdim/nsx;
	k_x = dvector(1,nns);
    k_y = dvector(1,nns);
    if (old==1){
    	FToldRE= dvector(1,nns);
    	FToldIM= dvector(1,nns);
	}
    if (rough){
  	 	nx = dvector(1,nns);
  	 	ny = dvector(1,nns);
		dstrike=dvector(1,nns);
		ddip=dvector(1,nns);
		drake=dvector(1,nns);
		REsurf   = dvector(1,nns);
		IMsurf   = dvector(1,nns);
		FTsurfRE = dvector(1,nns);
		FTsurfIM = dvector(1,nns);
    }
	REold    = dvector(1,nns);
	rakes    = dvector(1,nns);
	IMold    = dvector(1,nns);
	dipold   = dvector(1,nns);
	rakeold   = dvector(1,nns);
	strikeold   = dvector(1,nns);
	REslip   = dvector(1,nns);
	IMslip   = dvector(1,nns);
	slip   	= dvector(1,nns);		//final slip(after tapering, rescaling).
	FTrandRE = dvector(1,nns);
	FTrandIM = dvector(1,nns);
	FTrand2RE = dvector(1,nns);
	FTrand2IM = dvector(1,nns);
	specRE   = dvector(1,nns);
	specIM   = dvector(1,nns);

	minslip0=1e30, maxslip0=-1e30;
	for (int np=1; np<=ns; np++){
		if (old){
			REold[np]=pow((*eqkfm1).slip_str[np]*(*eqkfm1).slip_str[np]+(*eqkfm1).slip_dip[np]*(*eqkfm1).slip_dip[np],0.5);
			rakes[np]=(-180/pi)*atan2((*eqkfm1).slip_dip[np],(*eqkfm1).slip_str[np]);
			minslip0= fmin(minslip0, REold[np]);
			maxslip0= fmax(maxslip0, REold[np]);
		}
		else REold[np]=0.0;
		IMold[np]=0;
	}

	if ((old) & ((maxslip0-minslip0)/fmax(fabs(maxslip0),fabs(minslip0))<0.01))	{
		print_screen("Input slip model uniform in suomod1_hf --> model will be tapered before Fourier Transform.\n");
		print_logfile("Input slip model uniform in suomod1_hf --> model will be tapered before Fourier Transform.\n");
		if ((*eqkfm1).taper) suomod1_taper((*eqkfm1), eqkfm1);
		for (int np=1; np<=ns; np++) REold[np]=pow((*eqkfm1).slip_str[np]*(*eqkfm1).slip_str[np]+(*eqkfm1).slip_dip[np]*(*eqkfm1).slip_dip[np],0.5);
	}

	if (rough){
		for (int np=1;np<=ns;np++){
			if (ox[np]>maxox) maxox=ox[np];
			if (oy[np]>maxoy) maxoy=oy[np];
			if (ox[np]<minox) minox=ox[np];
			if (oy[np]<minoy) minoy=oy[np];
			rakeold[np]=-(180/pi)*atan2((*eqkfm1).slip_dip[np],(*eqkfm1).slip_str[np]);	//check this! (sign)
			if (isnan(rakeold[np])==1) rakeold[np]=rake;	//for example if slip=0 in both directions.
			strikeold[np]=strike;
			dipold[np]= dip1;
		}

		start_x=minox-0.5*odiscx+0.5*disc;
		start_y=minoy-0.5*odiscy+0.5*disc;

		maxoy=maxoy+odiscy;
		i=1;
		for (k=1; k<=nsy; k++){
			for (j=1;j<=nsx;j++){
				nx[i]=start_x+(j-1)*disc; //relative position along strike
				ny[i]=start_y+(k-1)*disc; //relative position along dip
				i++;
			}
		}
	}

    // generate wavenumber vector

	kny_x=1./disc/2.;  // Nyquist wave numbers=maximal wave number
	kny_y=1./disc/2.;
	dk_x=1./oxdim;  // discretisation in wavenumber space
	dk_y=1./oydim;
	for (int i=0;i<nns;i++){
		 l=floor(i/nsx);
		 if (l>nsy/2) l=floor(i/nsx)-nsy;
		 j=i%nsx;
		 if (j>nsx/2) j=i%nsx-nsx;
		 k_x[i+1]=j*dk_x;
		 k_y[i+1]=l*dk_y;
	}
    
    // 2dft of old distribution if given
	if (old) dft2d(nsy,nsx,0, REold, IMold, FToldRE, FToldIM);
    
    //  create random field of size nns and make 2D-DFT of it
	for (i=1;i<=nns;i++){
		  REslip[i]=ran1(&seed2);
		  seed2=-seed2;
		  IMslip[i]=0;
	}
	dft2d(nsy,nsx,0, REslip, IMslip, FTrandRE, FTrandIM);

	if (rough){
	// create a second random field for the variation of strike and dip
		for (i=1;i<=nns;i++){
			REslip[i]=ran1(&seed2);
			seed2=-seed2;
			IMslip[i]=0;
		}
		dft2d(nsy,nsx,0, REslip, IMslip, FTrand2RE, FTrand2IM);
	}

	if (print_Fourier==1){
		if (old){
			sprintf(fname,"%s/Fourier_old.out", logfolder);
			fout=fopen(fname,"w");
			if (fout==NULL) print_screen("Warning: could not write output file for suomod1_hf.\n");
			else{
				for (k=1;k<=nns;k++) fprintf(fout,"%lf\t %lf\t%lf\t%lf\n", k_x[k], k_y[k], FToldRE[k], FToldIM[k]);
				fclose(fout);
			}
		}

		sprintf(fname,"%s/Fourier_rand.out", logfolder);
		fout=fopen(fname,"w");
		if (fout==NULL) print_screen("Warning: could not write output file for suomod1_hf.\n");
		else{
			for (k=1;k<=nns;k++) fprintf(fout,"%lf\t %lf\t%lf\t%lf\n", k_x[k], k_y[k], FTrandRE[k], FTrandIM[k]);
			fclose(fout);
		}
	}
    
    //  calculate spectrum with random phases and amplitude decay 1/(k)^(1+H) corresponding to Hurst exponent given above for k higher than the corner wavenumber kc

	kc=pow(10,(1.82-0.5*(Me))); 		// corner wavenumber according to Causse et al. (2010),

	kcx=2.5/oxdim;
	kcy=2.5/oydim;
	for (i=1;i<=nns;i++){
		k_i=sqrt(pow(fabs(k_x[i]), 2)+ pow(fabs(k_y[i]), 2));
//		A= (k_i< kc) ? pow(fabs(kc), H+1.0) : pow(fabs(k_x[i]), H+1.0)+ pow(fabs(k_y[i]), H+1.0);
//		specRE[i]=FTrandRE[i]/A;
//		specIM[i]=FTrandIM[i]/A;

		A= (k_i< kc) ? fabs(kc) : fabs(k_i);
		specRE[i]=FTrandRE[i]*pow(A,-(1.0+H));
		specIM[i]=FTrandIM[i]*pow(A,-(1.0+H));

		if (rough){
			A= (k_i< kc) ? pow(fabs(kc), 2.4) : pow(fabs(k_x[i]), H+2.4)+ pow(fabs(k_y[i]), 2.4);
			FTsurfRE[i]=FTrand2RE[i]/A;
			FTsurfIM[i]=FTrand2IM[i]/A;
		}
	}


	if (print_Fourier==1){
		sprintf(fname,"%s/Fourier_filt.out", logfolder);
		fout=fopen(fname,"w");
		if (fout==NULL) print_screen("Warning: could not write output file for suomod1_hf.\n");
		else{
			for (k=1;k<=nns;k++) fprintf(fout,"%lf\t %lf\t%lf\t%lf\n", k_x[k], k_y[k],specRE[k], specRE[k]);
			fclose(fout);
		}
	}

	// combine spectra of old and random distribution if old slip is given
	if (old==1){
		l=0;
		diffspec=0;
		for(i=1;i<=nns;i++){
			if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kcy){
			  diffspec+= log10(sqrt(specRE[i]*specRE[i]+specIM[i]*specIM[i]))-log10(sqrt(FToldRE[i]*FToldRE[i]+FToldIM[i]*FToldIM[i]));
			  l++;
			}
		}
		diffspec=diffspec/l;
		for(i=1;i<=nns;i++){
			if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kcy){
				specRE[i]=FToldRE[i]*pow(10,diffspec);
				specIM[i]=FToldIM[i]*pow(10,diffspec);
			}
		}
	}


	// inverse fourier transform to obtain final slip distribution:
	dft2d(nsy,nsx,1, specRE, specIM, REslip, IMslip);
	// set negative values of slip to 0 by substracting the minimum value:
    minslip=REslip[1];
    for (i=1;i<=nns;i++) if (REslip[i] < minslip) minslip=REslip[i];

    if (rough) {
		dft2d(nsy,nsx,1, FTsurfRE, FTsurfIM, REsurf, IMsurf);
	    minsurf=REsurf[1];
	    for (i=1;i<=nns;i++) if (REsurf[i] < minsurf) minsurf=REsurf[i];

	    // adjust height of fractal surface to aspect ratio
	    maxsurf=REsurf[1]-minsurf;
		for (i=1;i<=nns;i++)
		{ REsurf[i]=REsurf[i]-minsurf;
		  if (REsurf[i]>maxsurf)
			{maxsurf=REsurf[i];
			}
		}
		if (oxdim > oydim)
		{ height= oydim*aspect;
		}
		else	{ height= oxdim*aspect;
		}
		for (i=1;i<=nns;i++)
			  { REsurf[i]=REsurf[i]/maxsurf*height;
			  }
	  // calculate deviations of strike and dip as local angles along the fractal surface
		for (i=1;i<=nns;i++)
		  { if (nx[i]==0 )
				  { dstrike[i]=strikeold[i]+atan((REsurf[i+1]-REsurf[i])/disc)*180.0/pi;
				  }
			else if (nx[i]==(nsx-1)*disc)
				  { dstrike[i]=strikeold[i]+atan((REsurf[i]-REsurf[i-1])/disc)*180.0/pi;
				  }
			else  { dstrike[i]=strikeold[i]+atan((REsurf[i+1]-REsurf[i-1])/2/disc)*180.0/pi;
				  }
			if (ny[i]==0)
				  { ddip[i]=dipold[i]+atan((REsurf[i+(nsx)]-REsurf[i])/disc)*180.0/pi;
				  }
			else if (ny[i]==(nsy-1)*disc)
				  { ddip[i]=dipold[i]+atan((REsurf[i]-REsurf[i-nsx])/disc)*180.0/pi;
				  }
			else  { ddip[i]=dipold[i]+atan((REsurf[i+nsx]-REsurf[i-nsx])/2/disc)*180.0/pi;
				  }
			drake[i]=rakeold[i]+REsurf[i]/height*20.0-10.0;	//degrees or radians? what is this formula?
		  }
	}


   /* old stuff, not needed here.

    /// determine a possible hypocentre: find patch with maximum slip, determine the nearest point at which slip = 1/2max

    maxslip=(REslip[1]-(minslip))*taper[1]+0.001;
      for (i=1;i<=nns;i++)
    	{ if (maxslip<(REslip[i]-(minslip))*taper[i]+0.001)
    		{ maxslip=(REslip[i]-(minslip))*taper[i]+0.001;
    		  maxlat=nx[i];
    		  maxlon=ny[i];
    		}
    	}
      distance=oxdim;
      for (i=1;i<=nns;i++){
          dist=sqrt((nx[i]-maxlat)*(nx[i]-maxlat)+(ny[i]-maxlon)*(ny[i]-maxlon));
    	  if ((0.5*maxslip > ((REslip[i]-(minslip))*taper[i]+0.001)) && (distance>dist) && taper[i]>0.95){
    		  lat1=lat[i];
    		  lon1=lon[i];
    		  xhypo=nx[i];
    		  yhypo=ny[i];
    		  distance=dist;
    	  }
      }
    //// Calculate variation for rupture velocity
      vel=dvector(1,nns);
      for (i=1;i<=nns;i++) vel[i]=((REslip[i]-(minslip))*taper[i]+0.001)/maxslip*velmean*fraction*2+velmean-velmean*fraction;

     */


    //-------------------------------------Fill in eqkfm2.---------------------------------------//

	copy_eqkfm_noslipmodel((*eqkfm1), eqkfm2);

	(*eqkfm2).L=oxdim;
	(*eqkfm2).W=oydim;
	(*eqkfm2).np_st=nsx;
	(*eqkfm2).np_di=nsy;
	(*eqkfm2).pos_s=(*eqkfm1).pos_s;
	(*eqkfm2).pos_d=(*eqkfm1).pos_d;
	if ((*eqkfm2).slip_str==NULL) (*eqkfm2).slip_str=dvector(1,nns);
	if ((*eqkfm2).slip_dip==NULL) (*eqkfm2).slip_dip=dvector(1,nns);

 	(*eqkfm2).noise=1;



 	//----------set up vectors containing slip-----------//

	for (int p=1; p<=nns; p++) slip[p]=REslip[p];		//final slip.
	scale_to_mag((*eqkfm1), eqkfm2, slip,rakes);

	if (noise_only) {
		for (int p=1; p<=nns; p++){
			(*eqkfm2).slip_dip[p]-=(*eqkfm1).slip_dip[p];
			(*eqkfm2).slip_str[p]-=(*eqkfm1).slip_str[p];
		}
	}

	if ((*eqkfm2).taper) suomod1_taper(*eqkfm2, eqkfm2);

	if (noise_only) {
		(*eqkfm2).tot_slip=(*eqkfm1).tot_slip;
	}

 	//-------------print out--------------//


	if (print_Fourier==1){
		for (k=1;k<=nns;k++) {
			FTrandRE[k]=0.0;
			FTrandIM[k]=0.0;
		}
		dft2d(nsy,nsx,0, slip, IMold, FTrandRE, FTrandIM);
		sprintf(fname,"%s/Fourier.out", logfolder);
		fout=fopen(fname,"w");
		if (fout==NULL) print_screen("Warning: could not write output file for suomod1_hf.\n");
		else{
			for (k=1;k<=nns;k++) fprintf(fout,"%lf\t %lf\t %lf\t %lf\t\n", k_x[k],k_y[k], FTrandRE[k], FTrandIM[k]);
			fclose(fout);
		}
	}

	if (printout==1) {
		sprintf(outname2, "%s/source_hf.out",logfolder);
		print_slipmodel(outname2,eqkfm2,1);
	}

  //-------------free memory--------------//

	if (old==1){
		free_dvector(FToldRE,1,nns);
		free_dvector(FToldIM,1,nns);
	}
	free_dvector(REold,1,nns);
	free_dvector(IMold,1,nns);
	free_dvector(slip,1,nns);
	free_dvector(rakes,1,nns);
	free_dvector(k_x,1,nns);
	free_dvector(k_y,1,nns);
//	free_dvector(vel,1,nns);
	free_dvector(REslip,1,nns);
	free_dvector(IMslip,1,nns);
	free_dvector(FTrand2RE,1,nns);
	free_dvector(FTrand2IM,1,nns);
	free_dvector(FTrandRE,1,nns);
	free_dvector(FTrandIM,1,nns);
	free_dvector(specRE,1,nns);
	free_dvector(specIM,1,nns);
	free_dvector(strikeold,1,nns);
	free_dvector(dipold,1,nns);
	free_dvector(rakeold,1,nns);
	if (rough){
		free_dvector(dstrike,1,nns);
		free_dvector(ddip,1,nns);
		free_dvector(drake,1,nns);
		free_dvector(REsurf,1,nns);
		free_dvector(IMsurf,1,nns);
		free_dvector(FTsurfRE,1,nns);
		free_dvector(FTsurfIM,1,nns);
		free_dvector(nx,1,nns);
		free_dvector(ny,1,nns);
	}

	*seed = seed2;
	return(0);
}

int suomod1_addnoise(struct eqkfm eqkfm1, struct eqkfm eqkfm2, struct eqkfm *eqkfm){
	// sf is scale factor by which second slip model is scaled.
	//Normalization keeps magnitude of eqkfm1.

	int NP=eqkfm1.np_st*eqkfm1.np_di;
	double sf;	//scaling factor for noise.

	(*eqkfm).lat=eqkfm1.lat;
	(*eqkfm).lon=eqkfm1.lon;
	(*eqkfm).depth=eqkfm1.depth;
	(*eqkfm).L=eqkfm1.L;
	(*eqkfm).W=eqkfm1.W;
	(*eqkfm).str1=eqkfm1.str1;
	(*eqkfm).dip1=eqkfm1.dip1;
	(*eqkfm).t=eqkfm1.t;
	(*eqkfm).str2=eqkfm1.str1;
	(*eqkfm).dip2=eqkfm1.dip1;
	(*eqkfm).whichfm=1;
	(*eqkfm).nsel=eqkfm1.nsel;
	(*eqkfm).selpoints= eqkfm1.selpoints;
	(*eqkfm).taper=eqkfm1.taper;
	(*eqkfm).np_st=eqkfm1.np_st;
	(*eqkfm).np_di=eqkfm1.np_di;
	(*eqkfm).pos_s=eqkfm1.pos_s;
	(*eqkfm).pos_d=eqkfm1.pos_d;
	(*eqkfm).mag=eqkfm1.mag;
	(*eqkfm).rake1=eqkfm1.rake1;
	(*eqkfm).rake2=eqkfm1.rake2;
    (*eqkfm).index_cat=eqkfm1.index_cat;	//index of event in catalog (only for the catalog used for LL calculation).
    (*eqkfm).slip_str=dvector(1,NP);
   	(*eqkfm).slip_dip=dvector(1,NP);
	(*eqkfm).is_slipmodel= eqkfm1.is_slipmodel;

   	double M0old, M0;

   	//find ratio of moments:
   	sf= eqkfm1.tot_slip/eqkfm2.tot_slip;

   	//normalize to keep moment constant:
   	M0old=M0=0;
   	for (int p=1; p<=NP; p++){
		(*eqkfm).slip_str[p]=eqkfm1.slip_str[p]+sf*eqkfm2.slip_str[p];
		(*eqkfm).slip_dip[p]=eqkfm1.slip_dip[p]+sf*eqkfm2.slip_dip[p];
		M0old+=sqrt(pow(eqkfm1.slip_dip[p],2)+pow(eqkfm1.slip_str[p],2));
		M0+=sqrt(pow((*eqkfm).slip_dip[p],2)+pow((*eqkfm).slip_str[p],2));
   	}
   	for (int p=1; p<=NP; p++){
   		(*eqkfm).slip_str[p]*=M0old/M0;
   		(*eqkfm).slip_dip[p]*=M0old/M0;
   	}

   	(*eqkfm).tot_slip=tot_slip(*eqkfm);
   	suomod1_taper((*eqkfm), eqkfm);
   	return(0);
}

int suomod1_cleanup(struct eqkfm *eqkfm2){

	int nns=(*eqkfm2).np_di*(*eqkfm2).np_st;

	free_dvector((*eqkfm2).pos_s,1,nns);
	free_dvector((*eqkfm2).pos_d,1,nns);
	free_dvector((*eqkfm2).slip_str,1,nns);
	free_dvector((*eqkfm2).slip_dip,1,nns);

	return 0;
}
