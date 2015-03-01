/*
 * grid_variability.c
 *
 *  Created on: Nov 4, 2013
 *      Author: camcat
 */

#include "grid_variability.h"

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../geom/coord_trafos.h"
#include "../inp_out/print_output.h"
#include "../inp_out/read_crust.h"
#include "../inp_out/read_eqkfm.h"
#include "../okada/okadaDCFS.h"
#include "../seis/soumod1.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"

int grid_variability(){

/* Compares the variaiblity calculated using simplified method (i.e. comparing each cell center with its neighbouring cells)
 * to the variability estimate more precisely, i.e. by sampling a large number (Nsub) of random locations in the cell.
 * Loops over a set of simple slip models with large variabiltiy (H=5?).
 */

	double H=0.1;	//Hurst coefficient.
	char *testfolder="test/grid_var";
	char *crust_file="input/Tohoku_simple_vert.inp";
	int N=2;
	int res_ratio=5;
	long seed=-19329935;
	double res=1.5;
	struct eqkfm eqfm;
	struct crust crst, crst2;
	FILE *fout0, *fout;
	char fname[120], fname0[120];
	struct pscmp dcfs;
	double *cmb0, **interp_DCFS;
	double lat0=-38.5, lon0=142.0, dep0=0.0;
	int Nlat=16, Nlon=16, Ndep=4, NP=Nlat*Nlon*Ndep, NP2;
	//int Nlat=5, Nlon=5, Ndep=1, NP=Nlat*Nlon*Ndep, NP2;
	int Nlat2=Nlat*res_ratio, Nlon2=Nlon*res_ratio, Ndep2=Ndep*res_ratio;
	double *latgrid, *longrid, *depgrid;
	double *latgrid2, *longrid2, *depgrid2;
	int *allp, *allp2;
	double dlat=0.8, dlon=0.8, ddep=4.0;
	int P[3], Pn[3], nsub[3];
	int Nlatlon=Nlat*Nlon;	//no of points per dimension (old geometry);
	int ii;
	int **map;

	double strikes[]={90, 90};
	double dips[]={90, 15};
	double rakes[]={0, 90};


	//low res grid:
	latgrid=dvector(1,NP);
	longrid=dvector(1,NP);
	depgrid=dvector(1,NP);
	allp=ivector(1,NP);
	crst.x=dvector(1,NP);
	crst.y=dvector(1,NP);
	interp_DCFS=dmatrix(1,NP, 1,2);

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
	//read_pscmp_crust(crust_file,&crst);
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


	//high res grid:
	NP2=NP*pow(res_ratio,3);
	latgrid2=dvector(1,NP2);
	longrid2=dvector(1,NP2);
	depgrid2=dvector(1,NP2);
	allp2=ivector(1,NP2);
	crst2.x=dvector(1,NP2);
	crst2.y=dvector(1,NP2);

	for (int i=1; i<=NP2; i++) allp2[i]=i;

	for (int i=1; i<=Nlat2; i++){
		for (int j=1; j<=Nlon2; j++){
			for (int k=1; k<=Ndep2; k++){
				latgrid2[Nlat2*Nlon2*(k-1)+Nlat2*(j-1)+i]= lat0-0.5*dlat+(dlat/Nlat2)*(i-0.5);
				longrid2[Nlat2*Nlon2*(k-1)+Nlat2*(j-1)+i]= lon0-0.5*dlon+(dlon/Nlon2)*(j-0.5);
				depgrid2[Nlat2*Nlon2*(k-1)+Nlat2*(j-1)+i]= (ddep/Ndep2)*(k-0.5);
			}
		}
	}

	//broken: function does not exist anymore
	//read_pscmp_crust(crust_file,&crst2);
	crst2.nLat=Nlat2;
	crst2.nLon=Nlon2;
	crst2.nD=Ndep2;
	crst2.lat0=lat0;
	crst2.lon0=lon0;
	crst2.lat=latgrid2;
	crst2.lon=longrid2;
	crst2.depth=depgrid2;
	crst2.dlat=dlat/Nlat2;
	crst2.dlon=dlon/Nlon2;
	crst2.ddepth=ddep/Ndep2;
	for (int i=1; i<=NP2; i++) latlon2localcartesian(latgrid2[i], longrid2[i], lat0, lon0, crst2.y+i, crst2.x+i);

	//---------map between low, high res:

	map=imatrix(1,NP, 1, pow(res_ratio,3));

	//no. of grid points between high resolution geometry and low resolution geometry (per dimension);
	nsub[0]=res_ratio;
	nsub[1]=res_ratio;
	nsub[2]=res_ratio;

	for (int i=1; i<=NP; i++){
		//reshape linear array into 3x3 array: i -> (P1,P2,P3).
		ii=0;
		P[0]=((i-1)%(Nlatlon))%Nlat+1;
		P[1]=(((i-1)%Nlatlon)+1-P[0])/Nlat +1;
		P[2]=(i-P[0]-(P[1]-1)*Nlat)/Nlatlon +1;

		for (int x=1; x<=nsub[0]; x++){
			Pn[0]=(P[0]-1)*nsub[0]+x;
			for (int y=1; y<=nsub[1]; y++){
				Pn[1]=(P[1]-1)*nsub[1]+y;
				for (int z=1; z<=nsub[2]; z++){
					ii+=1;
					Pn[2]=(P[2]-1)*nsub[2]+z;
					map[i][ii]=Pn[0]+Nlat2*(Pn[1]-1)+Nlat2*Nlon2*(Pn[2]-1);
				}
			}
		}
	}


	////////////////////////////////////////////////////////////////////

	sprintf(fname0, "%s/fm_slipmodels_summary.txt",testfolder);
	fout0=fopen(fname0,"w");

	fprintf(fout0, "strike\t dip \t rake \t magnitude \t resolution \n");
	eqfm.lat=lat0;
	eqfm.lon=lon0;
	eqfm.whichfm=1;

	dcfs.NF=1;
	latlon2localcartesian(eqfm.lat, eqfm.lon, crst.lat0, crst.lon0, &(eqfm.y), &(eqfm.x));

	for (int n=0; n<N; n++){
		printf("n=%d\n",n);
		eqfm.depth=dep0;
		eqfm.str1=strikes[n];
		eqfm.dip1=dips[n];
		eqfm.rake1=rakes[n];
		eqfm.whichfm=1;
		eqfm.mag=6.0;
		//res=0.5+ran1(&seed);

		fprintf(fout0, "%.0lf\t%.0lf\t%.0lf\t%.1lf\t%.1lf\n", eqfm.str1, eqfm.dip1, eqfm.rake1, eqfm.mag, res);
		focmec2slipmodel(crst, &eqfm, res, 1, 1);
		//suomod1_hf(eqfm, &eqfm, H, &seed, 0);	//fixme broke: function does not exist anymore

		sprintf(fname, "%s/slipmodel_%d.dat",testfolder, n);
		print_slipmodel(fname, &eqfm, 1);

		dcfs.nsel=eqfm.nsel=NP;
		dcfs.which_pts=eqfm.selpoints=allp;
		dcfs.S=d3tensor(1,dcfs.nsel, 1,3,1,3);
		dcfs.cmb=dvector(1,dcfs.nsel);

//		okadaDCFS(dcfs, &eqfm, 1, crst, &eqfm.str1, &eqfm.dip1, 0);
		cmb0=dcfs.cmb;
		interp_nn(NP, Nlat, Nlon, Ndep,dcfs.cmb, interp_DCFS, 0, NULL);

		sprintf(fname, "%s/cmb_%d_lowres.dat",testfolder, n);
		print_grid(fname, dcfs, crst, NULL);

		dcfs.nsel=eqfm.nsel=NP2;
		dcfs.which_pts=eqfm.selpoints=allp2;
		dcfs.S=d3tensor(1,dcfs.nsel, 1,3,1,3);
		dcfs.cmb=dvector(1,dcfs.nsel);

//		okadaDCFS(dcfs, &eqfm, 1, crst2, &eqfm.str1, &eqfm.dip1, 0);
		sprintf(fname, "%s/cmb_%d_highres.dat",testfolder, n);
		print_grid(fname, dcfs, crst2, NULL);

		sprintf(fname, "%s/cmb_%d_info.dat",testfolder, n);
		fout=fopen(fname,"w");
		for (int i=1; i<=NP; i++){
			fprintf(fout, "%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t", crst.lat[i], crst.lon[i], crst.depth[i], cmb0[i], interp_DCFS[i][1], interp_DCFS[i][2]);
			fprintf(fout,"\n");
		}
		fclose(fout);

		sprintf(fname, "%s/cmb_%d_cmb.dat",testfolder, n);
		fout=fopen(fname,"w");
		for (int i=1; i<=NP; i++){
			for (int k=1; k<=pow(res_ratio,3); k++) fprintf(fout, "%.3lf\t", dcfs.cmb[map[i][k]]);
			fprintf(fout,"\n");
		}
		fclose(fout);

		sprintf(fname, "%s/cmb_%d_lats.dat",testfolder, n);
		fout=fopen(fname,"w");
		for (int i=1; i<=NP; i++){
			for (int k=1; k<=pow(res_ratio,3); k++) fprintf(fout, "%.3lf\t", latgrid2[map[i][k]]);
			fprintf(fout,"\n");
		}
		fclose(fout);

		sprintf(fname, "%s/cmb_%d_lons.dat",testfolder, n);
		fout=fopen(fname,"w");
		for (int i=1; i<=NP; i++){
			for (int k=1; k<=pow(res_ratio,3); k++) fprintf(fout, "%.3lf\t", longrid2[map[i][k]]);
			fprintf(fout,"\n");
		}
		fclose(fout);

		sprintf(fname, "%s/cmb_%d_deps.dat",testfolder, n);
		fout=fopen(fname,"w");
		for (int i=1; i<=NP; i++){
			for (int k=1; k<=pow(res_ratio,3); k++) fprintf(fout, "%.3lf\t", depgrid2[map[i][k]]);
			fprintf(fout,"\n");
		}
		fclose(fout);
	}

	fclose(fout0);

	return 0;
}

