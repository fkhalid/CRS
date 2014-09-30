/*
 * struct_conversions.c
 *
 *  Created on: Dec 22, 2011
 *      Author: camcat
 */

#include "struct_conversions.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

//------------------ combining -------------------//

int *combine_eqkfm(struct eqkfm *eqkfm1, struct eqkfm *eqkfm2, int N1, int N2,
				   double tend, double dt, double dM, double dR, int overwrite) {
//joins 2 eqkfm catalogs, where each member corresponds to an earthquake (i.e. no multifault events). Keeps event number of first catalog.
//dt, dM are ranges within which earthquakes are considered to be the same.
//indices range between [0,N1-1] and [0,N2-1].

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int n1=0, n10=0, n12=0; //indices of next and previous event and closest in time.
	double dist20, dist2, dlon;
	int selected, *sel, *sel1;
	int N20=N2, N10=N1;
	double dx, dy, dz, r;
	int outfile=0, not_selected=0;	//outfile can be reactivated for printing out selected events (e.g. for debugging).
	char fname[120];
	FILE *fout;

	//recalculate N1, N2 to exclude events after tend.
	N2=N1=0;
	while (N2<N20 && eqkfm2[N2].t<=tend) N2++;
	while (N1<N10 && eqkfm1[N1].t<=tend) N1++;

	if (N1==0 | N2==0) {
		print_screen("Warning - one of the eqkfm structures is empty! (combine_eqkfm) \n");
		print_logfile("Warning - one of the eqkfm structures is empty! (combine_eqkfm) \n");
		return NULL;
	}

	sel=ivector(0,N2-1);
	sel1=ivector(0,N1-1);
	for (int n=0; n<N1; n++) sel1[n]=0;

	for (int n2=0; n2<N2; n2++){
		selected=0;

		n10=n12;
		while (n10<N1-1 && eqkfm1[n10].t<=eqkfm2[n2].t-dt) n10++;
		n1=n10;
		while (n1<N1-1 && eqkfm1[n1].t<=eqkfm2[n2].t+dt) n1++;

		dist2=1e30;
		for (int n=n10; n<=n1; n++){
			if (sel1[n]!=0) continue;	//todo should not assign elements on a first come, first served basis...
			dist20=pow((eqkfm2[n2].t-eqkfm1[n].t)/dt,2)+pow((eqkfm2[n2].mag-eqkfm1[n].mag)/dM,2);
			if (dist2>dist20){
				dist2=dist20;
				n12=n;
			}
		}

		dx=Re*(eqkfm2[n2].lat-eqkfm1[n12].lat)*DEG2RAD;
		dlon=eqkfm2[n2].lon-eqkfm1[n12].lon;
		if (fabs(dlon)>180) dlon=(dlon>0) ? dlon-360 : dlon+360;
		dy=Re*(eqkfm2[n2].lon-eqkfm1[n12].lon)*DEG2RAD*cos(eqkfm1[n12].lat*DEG2RAD);
		dz=eqkfm2[n2].depth-eqkfm1[n12].depth;
		r=sqrt(dx*dx+dy*dy);

		if (fabs(eqkfm2[n2].mag-eqkfm1[n12].mag) <=(dM+0.001) && fabs(eqkfm2[n2].t-eqkfm1[n12].t) <=dt && r<=dR){
			if(procId == 0) {
				if (outfile) fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",eqkfm2[n2].t, eqkfm2[n2].lat, eqkfm2[n2].lon, eqkfm2[n2].depth, eqkfm2[n2].mag, 1, eqkfm1[n12].t, eqkfm1[n12].lat, eqkfm1[n12].lon, eqkfm1[n12].depth, eqkfm1[n12].mag,r);
			}
			if (overwrite==1){
				//copy_eqkfm_nolocation_noindex_notime(eqkfm2[n2], eqkfm1+n12);	//todo decide which one is better (temporary change for iquique).
				copy_eqkfm_noindex_notime(eqkfm2[n2], eqkfm1+n12);
			}
			selected+=1;
			sel[n2]=n12;
			sel1[n12]=1;
		}

		else {
			if(procId == 0) {
				if (outfile) fprintf(fout,"%.8lf\t%lf\t%lf\t%lf\t%lf\t%d\t%.8lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",eqkfm2[n2].t, eqkfm2[n2].lat, eqkfm2[n2].lon, eqkfm2[n2].depth, eqkfm2[n2].mag, 0, 0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0);
			}
		}

		if (selected!=1 && eqkfm2[n2].t<tend) {
			sel[n2]=-1;
			if (!selected) {
				not_selected+=1;
				if (extra_verbose) {
					print_screen("Warning: element %d [t=%lf, Mw=%lf, d=%.3lf] from eqkfm2 missing in eqkfm1 (function: combine_eqkfm)!!\n",n2,eqkfm2[n2].t,eqkfm2[n2].mag, eqkfm2[n2].depth);
					print_logfile("Warning: element %d [t=%lf, Mw=%lf, d=%.3lf] from eqkfm2 missing in eqkfm1 (function: combine_eqkfm)!!\n",n2,eqkfm2[n2].t,eqkfm2[n2].mag, eqkfm2[n2].depth);
				}
			}
		}
	}
			
	if (not_selected){
		print_screen("Warning: %d elements from focal mechanism catalog missing from earthquake catalog.\n", not_selected);
		print_logfile("Warning: %d elements from focal mechanism catalog missing from earthquake catalog.\n", not_selected);
	}

	if(procId == 0) {
		if (outfile) fclose(fout);
	}

	return sel;
}

int *combine_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM){
	//t1, mX have indices [0...NX-1].
	int n1=0, n10=0, n12=0; //indices of next and previous event and closest in time.
	int selected;
	int *sel, *sel1;
	double dist2, dist20;

	if (!N2) return NULL;
	if (!N1) {
		sel= ivector(0,N2-1);
		for (int n=0; n<N2; n++) sel[n]=-1;
		return sel;
	}

	sel= ivector(0,N2-1);
	sel1= ivector(0,N1-1);
	for (int n=0; n<N1; n++) sel1[n]=0;

	for (int n2=0; n2<N2; n2++){
		selected=0;

		n10=n12;
		while (n10<N1-1 && t1[n10]<=t2[n2]-dt) n10++;
		n1=n10;
		while (n1<N1-1 && t1[n1]<=t2[n2]+dt) n1++;

		dist2=1e30;
		for (int n=n10; n<=n1; n++){
			if (sel1[n]!=0) continue;
			dist20=	pow((t2[n2]-t1[n])/dt,2)+pow((m2[n2]-m1[n])/dM,2);
			if (dist2>dist20){
				dist2=dist20;
				n12=n;
			}
		}

		if (fabs(m2[n2]-m1[n12]) <=(dM+0.001) && fabs(t2[n2]-t1[n12]) <=dt){
			selected+=1;
			sel[n2]=n12;
			sel1[n12]=1;
		}
		else sel[n2]=-1;

	}
	return sel;
}

double **union_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM, int ***ind, int *tot){
	//t1, mX have indices [0...NX-1].
	//gives times and magnitude combined (also non common elements).
	//value of -1 in ind[x][y] means that element y was not found in one [tx mx], otherwise index is given.
	//results and ind have indices [0...tot-1].

	int n1=0, n10=0, n12=0, n120; //indices of next and previous event, closest and closest to previous element.
	int selected, count=0;
	int *sel, *sel1;
	double dist2, dist20;
	double **res;

	if (ind) *ind=imatrix(1,2,0,N1+N2);
	res=dmatrix(1,2,0,N1+N2);

	if (!N1){
		if (tot) *tot=N2;
		for (int i=0; i<N2; i++){
			if (ind){
				(*ind)[1][i]=-1;
				(*ind)[2][i]=i;
			}
			res[1][i]=t2[i];
			res[2][i]=m2[i];
		}
		return res;
	}

	if (!N2){
		if (tot) *tot=N1;
		for (int i=0; i<N1; i++){
			if (ind){
				(*ind)[1][i]=i;
				(*ind)[2][i]=-1;
			}
			res[1][i]=t1[i];
			res[2][i]=m1[i];
		}
		return res;
	}


	sel=ivector(0,N2-1);
	sel1=ivector(0,N1-1);

	for (int n=0; n<N1; n++) sel1[n]=0;

	for (int n2=0; n2<N2; n2++){
		selected=0;

		n10=n120=n12;
		while (n10<N1-1 && t1[n10]<=t2[n2]-dt) n10++;
		n1=n10;
		while (n1<N1-1 && t1[n1]<=t2[n2]+dt) n1++;

		dist2=1e30;
		for (int n=n10; n<=n1; n++){
			if (sel1[n]!=0) continue;
			dist20=	pow((t2[n2]-t1[n])/dt,2)+pow((m2[n2]-m1[n])/dM,2);
			if (dist2>dist20){
				dist2=dist20;
				n12=n;
			}
		}

		if (fabs(m2[n2]-m1[n12]) <=(dM+0.001) && fabs(t2[n2]-t1[n12]) <=dt){
			selected+=1;
			sel[n2]=n12;
			sel1[n12]=1;
		}
		else sel[n2]=-1;

		int n0=n120;
		while (n0<n12 && t1[n0]<t2[n2]){
			if(sel1[n0]==0) {
				res[1][count]=t1[n0];
				res[2][count]=m1[n0];
				if (ind!= (int ***) 0){
					(*ind)[1][count]=n0;
					(*ind)[2][count]=-1;
				}
				count+=1;
			}
			n0+=1;
		}
		res[1][count]=t2[n2];
		res[2][count]=m2[n2];
		if (ind!= (int ***) 0){
			(*ind)[1][count]=sel[n2];
			(*ind)[2][count]=n2;
		}
		count+=1;
		for	(int n=n0; n<n12; n++){
			if (sel1[n]==0){
				res[1][count]=t1[n];
				res[2][count]=m1[n];
				if (ind!= (int ***) 0){
					(*ind)[1][count]=n;
					(*ind)[2][count]=-1;
				}
				count+=1;
			}
		}
	}

	if (sel1[n12]==1) n120=n12+1;
	else n120=n12;

	for (int n=n120; n<N1; n++){
		res[1][count]=t1[n];
		res[2][count]=m1[n];
		if (ind!= (int ***) 0){
			(*ind)[1][count]=n;
			(*ind)[2][count]=-1;
		}
		count+=1;
	}

	*tot=count;
	return res;
}


//------------------ filtering -------------------//

void eqk_filter(struct eqkfm **eqkfm1, int *Ntot, double Mag, double Depth){
	//inefficient (3 loops), but uses as little memory as possible.
	//if free==1, frees memory.
	struct eqkfm *eqkfm0;
	int j=0;
	int Ntot_new=0;

	for (int i=0; i<(*Ntot); i++){
		if ((*eqkfm1)[i].mag>=Mag && (*eqkfm1)[i].depth<=Depth) Ntot_new+=1;
	}

	eqkfm0 = eqkfm_array(0, Ntot_new-1);

	for (int i=0; i<(*Ntot); i++){
		if ((*eqkfm1)[i].mag>=Mag && (*eqkfm1)[i].depth<=Depth){
			copy_eqkfm_all((*eqkfm1)[i],eqkfm0+j);
			j+=1;
		}
	}

	*eqkfm1 = eqkfm_array(0, Ntot_new-1);

	for (int i=0; i<Ntot_new; i++) copy_eqkfm_all(eqkfm0[i],(*eqkfm1)+i);
	*Ntot=Ntot_new;
	print_screen("%d events with Mw>=%.3lf, z<=%.3lf selected from eqkfm (eqk_filter).\n",Ntot_new, Mag, Depth);
	print_logfile("%d events with Mw>=%.3lf, z<=%.3lf selected from eqkfm (eqk_filter).\n",Ntot_new, Mag, Depth);
	return;
}

//--------------------extracting 1d arrays------------------------//

double *timesfromeqkfm(struct eqkfm *eqkfm1, int N, int *NF){
/* simply copy times from eqkfm to double vector. indices: [0...N]
 * NF contains no. of faults for each event; if NULL, assume single fault events. */


	double *times=dvector(0,N-1);
	int counter=0;

	for (int i=0; i<N; i++) {
		times[i]=eqkfm1[counter].t;
		counter= (NF==NULL)? counter+1 : counter+NF[i];
	}
	return times;

}

double *magssfromeqkfm(struct eqkfm *eqkfm1, int N, int *NF){
/*simply copy times from eqkfm to double vector. indices: [0...N]
 * NF contains no. of faults for each event; if NULL, assume single fault events. */

	double *mags=dvector(0,N-1);
	int counter=0, NF_i;
	double M0;
	for (int i=0; i<N; i++) {

		M0=0.0;
		NF_i= (NF==NULL)? 1 : NF[i];
		for (int k=0; k<NF_i; k++)	M0+=pow(10,1.5*(eqkfm1[counter+k].mag+6.0));
		mags[i]=(2.0/3.0)*log10(M0)-6.0;
		counter= counter+NF_i;
	}
	return mags;

}

double *timesfrompscmp(struct pscmp *DCFS, int N){
//simply copy times from DCFS to double vector. indices: [0...N]

double *times=dvector(0,N-1);
for (int i=0; i<N; i++) times[i]=DCFS[i].t;
return times;

}

double *magsfrompscmp(struct pscmp *DCFS, int N){
//simply copy magnitudes from DCFS to double vector. indices: [0...N]

double *mags=dvector(0,N-1);
for (int i=0; i<N; i++) mags[i]=DCFS[i].m;
return mags;

}

//todo move somewhere else?
void eqkfm2dist(struct eqkfm *eqkfm1, double *lats, double *lons, double *depths, int N, int Ntot, int all){
// if flag all=0, only find distance if eqkfm[i].is_slipmodel=0.

double x,y, *xs, *ys, Depth;
double lat0, lon0;
int nsel, pt;

ys=dvector(1,N);
xs=dvector(1,N);
lat0=0.5*(lats[N]+lats[1]);
lon0=0.5*(lons[N]+lons[1]);
for (int k0=1; k0<=N;k0++) latlon2localcartesian(lats[k0], lons[k0], lat0, lon0, ys+k0, xs+k0);

for (int i=0; i<Ntot; i++){
	if (all==1 || eqkfm1[i].is_slipmodel==0){
		nsel=eqkfm1[i].nsel;
		if (nsel==0) continue;

		latlon2localcartesian(eqkfm1[i].lat, eqkfm1[i].lon, lat0, lon0, &y, &x);
		Depth=eqkfm1[i].depth;
		eqkfm1[i].distance=dvector(1,nsel);	//TODO deallocate.

		for (int p=1; p<=nsel; p++) {
			pt=eqkfm1[i].selpoints[p];
			eqkfm1[i].distance[p]= sqrt(pow(ys[pt]-y,2)+pow(xs[pt]-x,2)+pow(depths[pt]-Depth,2));
		}
	}
}
}
