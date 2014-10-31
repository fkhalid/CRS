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

// todo [coverage] this block is never tested
double tot_slip (struct eqkfm eqfm1){

	double slip=0.0;
	int NP=eqfm1.np_di*eqfm1.np_st;

	for (int p=1; p<=NP; p++) {
		slip+=pow(eqfm1.slip_str[p]*eqfm1.slip_str[p]+eqfm1.slip_dip[p]*eqfm1.slip_dip[p],0.5);
	}

	return slip/NP;
}

// todo [coverage] this block is never tested
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

	#pragma omp parallel for
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

int suomod1_taper(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, int top, int bottom, int right, int left){
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

    // can be set to 1 to print extra files.
	printout=0;

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

    ns=nns=eqkfm1.np_st*eqkfm1.np_di;
    if (eqkfm1.np_st==1) left=right=0;
    if (eqkfm1.np_di==1) top=bottom=0;
    if (ns==1) {
    	copy_eqkfm_all(eqkfm1, eqkfm2);	//bug fix: using copy_eqkfm_slipmodel did not copy nsel, so no points were selected.
    	if (extra_verbose) {
    		print_screen("** Warning: model has a single patch, will not be tapered.**\n");
    	   	print_logfile("** Warning: model has a single patch, will not be tapered.**\n");
    	}
    	return 0;
    }

	// todo [coverage] this block is never tested
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

	#pragma omp parallel for
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

	#pragma omp parallel for
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

	#pragma omp parallel for
	for (int p=1; p<=nns; p++) {
		slip[p]=REslipo[p]*taper[p];		//final slip.
	}

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
