/*
 * setup_time.c
 *
 *  Created on: May 29, 2015
 *      Author: camcat
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "../defines.h"
#include "../util/merge.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"
#include "../util/splines_eqkfm.h"
#include "find_timesteps.h"
#include "lin_interp_eqkfm.h"


int timesteps_omori(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2,
		double smallstepstime, double TAU, double dtau, double timeTAU){

/* Calculates time steps based on the information in eqk_aft.
 * To obtain time steps with an increasing spacing, a function of the form t_{i}=t_{i-1}+K(t+c-teq)^p is used: for p=1 and a logarithmic stressing history,
 * this gives equal stresses within each time step. (NB: t refers to the start of each aseismic event).
 *
 * Input:
 *  t0, t1: start and end time.
 *  eqk_aft: array containing all the events to be considered.	Range [0...NFtot-1], where NFtot=sum(Nfaults);
 *  Nfaults:  number of faults per event. Range [0...NA-1].
 *
 *  smallstepstime: initial period during which smaller time steps are used (uses 0.3K instead of K).
 *	 K: parameter describing how closely spaced time steps are: t_{i}=t_{i-1}+K(t+c-teq)^p
 *
 * Output:
 *  times2: time steps. Memory is allocated here and the arrays is populated. Range [0...*L].
 *  L: largest index in times2.
 */

	int offset, L0, Ltot;
	int nev, nfaults;
	double Teq, tend;
	int err=0;
	double c=0.001, p=0.6;
	double K=0.6;

	//Dry run: find total no. of time steps, without filling out times2 (not allocated yet):
	offset=1;
	nev=nfaults=0;

	while(nev<NA && (*eqk_aft)[nfaults].t<t1){
		Teq=(*eqk_aft)[nfaults].t;
		tend= (nev<NA-1) ? fmin((*eqk_aft)[nfaults+Nfaults[nev]].t,t1) : t1;

		err+=findtimestepsomori(Teq, Teq, fmin(smallstepstime+Teq,tend), 0.3*K, 0.6, 0.001, NULL, NULL, L);
		if (smallstepstime+Teq<tend) {
			err+=findtimestepsomori(Teq, smallstepstime+Teq, tend, K, 0.6, 0.001, NULL, NULL, &L0);
			Ltot=offset=*L+L0+offset;
		}
		else{
			Ltot=offset=*L+offset;
		}
		nfaults+=Nfaults[nev];
		nev++;
	}

	//allocate memory for times2:
	*times2=dvector(0,Ltot);

	//run again, but fill in times2:
	offset=1;
	(*times2)[0]=fmin(t0,(*eqk_aft)[0].t)-1e-6;

	nev=nfaults=0;
	while(nev<NA && (*eqk_aft)[nfaults].t<t1){
		Teq=(*eqk_aft)[nfaults].t;
		tend= (nev<NA-1) ? fmin((*eqk_aft)[nfaults+Nfaults[nev]].t,t1) : t1;

		err+=findtimestepsomori(Teq, Teq, fmin(smallstepstime+Teq,tend), 0.3*K, 0.6, 0.001, (*times2)+offset, NULL, L);
		if (smallstepstime+Teq<tend) {
			err+=findtimestepsomori(Teq, smallstepstime+Teq, tend, K, 0.6, 0.001, (*times2)+*L+offset, NULL, &L0);
			Ltot=offset=*L+L0+offset;
		}
		else{
			Ltot=offset=*L+offset;
		}
		nfaults+=Nfaults[nev];
		nev++;
	}

	*L=Ltot;

	return err;

}

int timesteps_lin(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults,
		int *L, double **times2, int ***allind){

/* Calculates time steps by combining the snapshots times of all aseismic events, giving vector with combined and sorted snapshot times.
 * The total number of time steps is the sum of the time steps given for each event, and stressing histories are calculated accordingly.
 *
 * Input:
 *  t0, t1: start and end time.
 *  eqk_aft: array containing all the events to be considered.	Range [0...NFtot-1], where NFtot=sum(Nfaults);
 *  Nfaults:  number of faults per event. Range [0...NA-1]. *
 *
 * Output:
 *  times2: time steps. Memory is allocated here and the arrays is populated. Range [0...*L].
 *  L: largest index in times2.
 *  indices[n] contains the list of elements which originally belonged each of the events in eqk_aft. Memory allocated if a pointer to a NULL array is passed.
 */

	int nfaults=0;
	double **allts;	//contains all time steps lists.
	double *times2temp=NULL;	//temporary list containing time steps from snapshots.
	int *lens;	//lengths of time steps lists.

	allts=(double **) malloc((size_t)(NA*sizeof(double*)));
	lens= ivector(0,NA-1);

	//populate array containing all time steps:
	nfaults=0;
	for (int nev=0; nev<NA; nev++){
		allts[nev] = (double *) malloc((size_t)((*eqk_aft)[nfaults].nosnap+1)*sizeof(double));	//one extra element which corresponds to Teq (starting time when stresses are 0).
		allts[nev][0]=(*eqk_aft)[nfaults].t;
		for (int t=0; t<(*eqk_aft)[nfaults].nosnap; t++) allts[nev][t+1]=(*eqk_aft)[nfaults].ts[t];
		lens[nev]=(*eqk_aft)[nfaults].nosnap+1;
		nfaults+=Nfaults[nev];
	}

	//merge time steps into single array and keep indices:
	merge_multiple(allts, lens, NA, &times2temp, L, allind);

	//two extra elements at the start/end:
	(*times2)=dvector(0,*L+2);
	(*times2)[0]=fmin(t0, times2temp[0])-1e-6;
	(*times2)[*L+1]=fmax(t1, times2temp[*L-1])+1e-6;
	copy_vector(times2temp-1, times2, *L);
	*L+=2;

	for (int nev=0; nev<NA; nev++) free(allts[nev]);
	free(allts);

	return 0;
}

int setup_afterslip_multi_linear(double t0, double t1, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2){

/* Combines all stressing histories for all elements in eqk_aft, and rewrites them referred to combined time steps.
 *
 * Input:
 *  t0, t1: start and end time.
 *  eqk_aft: array containing all the events to be considered.	Range [0...NFtot-1], where NFtot=sum(Nfaults);
 *  Nfaults:  number of faults per event. Range [0...NA-1]. *
 *
 *
 * Output:
 *  times2: time steps. Memory is allocated here and the arrays is populated. Range [0...*L].
 *  L: largest index in times2.
 *
 *  The total number of time steps is the sum of the time steps given for each event, and stressing histories are calculated accordingly.
 *
 */

	int nfaults=0;
	int **allind=NULL;
	struct eqkfm *eq_aft= *eqk_aft;
	int printout_history=1, Nas; 	//can set to 1 to check if splines are giving correct stressing history.

	FILE *fout;
	char fname[120];

	if (printout_history){
		nfaults=0;
		for (int nev=0; nev<NA; nev++){
			Nas=(*eqk_aft)[nfaults].nosnap;
			sprintf(fname,"linear_old%d.dat",nev);
			fout=fopen(fname,"w");
			for (int l=0; l<Nas; l++) {
				fprintf(fout,"%.5e\t",(*eqk_aft)[nfaults].ts[l]);
			}
			fprintf(fout,"\n");
				for (int f=nfaults; f<nfaults+Nfaults[nev]; f++) {
					for (int p=1; p<=(*eqk_aft)[f].np_di*(*eqk_aft)[f].np_st; p++) {
						for (int l=0; l<Nas; l++) {
							if ((*eqk_aft)[f].allslip_str) fprintf(fout,"%.5e\t",(*eqk_aft)[f].allslip_str[l][p]);
						}
						fprintf(fout,"\n");
					}
				}
			fclose(fout);
			nfaults+=Nfaults[nev];
		}
	}

	//calculate combined time steps:
	timesteps_lin(t0, t1, eqk_aft, NA, Nfaults, L, times2, &allind);

	//recalculate slip at all time steps:
	nfaults=0;
	for (int nev=0; nev<NA; nev++){
		//need to shift indices by 1 since there is an extra element at the start:
		Nas=eq_aft[0].nosnap;
		for (int j=0; j<=Nas; j++) allind[nev][j]+=1;
		lin_interp_eqkfm(&eq_aft, Nfaults[nev], *L, *times2, allind[nev]);
		for (int f=0; f<Nfaults[nev]; f++) eq_aft[f].tevol=NULL;
		nfaults+=Nfaults[nev];
		eq_aft+=Nfaults[nev];
	}
	eq_aft-=nfaults;	//shift it back.

	//fixme check: should calculate differences here?

	if (printout_history){
		nfaults=0;
		for (int nev=0; nev<NA; nev++){
			sprintf(fname,"linear%d.dat",nev);
			fout=fopen(fname,"w");
			for (int l=0; l<=*L-1; l++) {
				fprintf(fout,"%.5e\t",(*times2)[l]);
			}
			fprintf(fout,"\n");
			for (int f=nfaults; f<nfaults+Nfaults[nev]; f++) {
				for (int p=1; p<=(*eqk_aft)[f].np_di*(*eqk_aft)[f].np_st; p++) {
					for (int l=0; l<=*L-1; l++) {
						if ((*eqk_aft)[f].allslip_str) fprintf(fout,"%.5e\t",(*eqk_aft)[f].allslip_str[l][p]);
					}
					fprintf(fout,"\n");
				}
			}
			fclose(fout);
			nfaults+=Nfaults[nev];
		}
	}

	return(0);

}

int setup_afterslip_single_linear(double t0, double t1, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2) {

//if splines==1, eqk_aft is substituted with more densily discretized version (L steps instead of Nas).
//Cs, ts, coefficients of temporal evolution functions (See Savage Parkfield paper). Nfun: no. of such funtions.
// eq_aft has indices: [0...Nas*Nfaults-1].
// logarithmic: flag indicating if stressing history is logarithmic. If not, use linear.

	int err=0;
//	double now, prev, norm, curr;
	int **allind=NULL;
	struct eqkfm *eq_aft;

	double times1[2];	//times for each aseismic event
	double **temp_tevol;
	int Nas=(*eqk_aft)[0].nosnap;
	int nfaults;

	eq_aft= *eqk_aft;

	//calculate combined time steps:
	timesteps_lin(t0, t1, eqk_aft, NA, Nfaults, L, times2, &allind);

	temp_tevol=dmatrix(1,1,0,*L-1);


	// Temporal evolution of afterslip.//
	nfaults=0;

	for (int nev=0; nev<NA; nev++){
		times1[0]=(*eqk_aft)[nfaults].t;
		times1[1]=(*eqk_aft)[nfaults].ts[0];

		//find value of tevol:
		fit_lin(times1, times2, Nas+1, *L, allind[nev], 1, NULL, temp_tevol);

		//assign tevol arrays:
		(*eqk_aft)[nfaults].tevol=dvector(1,*L);	//shifted by 1 because of copy_vector
		copy_vector(temp_tevol[1], &((*eqk_aft)[nfaults].tevol), *L);
		(*eqk_aft)[nfaults].tevol-=1;	//so it starts from 0
		for (int f=1; f<Nfaults[nev]; f++) {
			(*eqk_aft)[nfaults+f].tevol=(*eqk_aft)[nfaults].tevol;
		}
		nfaults+=Nfaults[nev];
	}

	free_dmatrix(temp_tevol, 1, 1, 0, *L-1);
	return(err!=0);
}

int setup_afterslip_multi_log(double t0, double t1, double *Cs, double *ts,
						 int Nfun, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2,
						 long *seed) {

//if splines==1, eqk_aft is substituted with more densily discretized version (L steps instead of Nas).
//Cs, ts, coefficients of temporal evolution functions (See Savage Parkfield paper). Nfun: no. of such funtions.
// eq_aft has indices: [0...Nas*Nfaults-1].

	int err=0;
	double TAU=200000, dtau=7000, timeTAU=183;	//todo allow to set from outside.
	double M0,mu;
	double smallstepstime=12;
	double now, prev, norm, curr;
	double Tendaft;	//Time to which cumulative afterslip snapshot refers.
	struct eqkfm *eq_aft;

	double Teq;
	double *t_afterslip;
	int Nas=(*eqk_aft)[0].nosnap;
	int nfaults;

	int printout_splines=1;	//can set to 1 to check if splines are giving correct stressing history.
	FILE *fout;
	char fname[120];

	eq_aft= *eqk_aft;

	err=timesteps_omori(t0, t1, eqk_aft, NA, Nfaults, L, times2, smallstepstime, TAU, dtau, timeTAU);

	if (printout_splines){
		nfaults=0;
		for (int nev=0; nev<NA; nev++){
			Nas=(*eqk_aft)[nfaults].nosnap;
			sprintf(fname,"splines_old%d.dat",nev);
			fout=fopen(fname,"w");
			for (int l=0; l<Nas; l++) {
				fprintf(fout,"%.5e\t",(*eqk_aft)[nfaults].ts[l]);
			}
			fprintf(fout,"\n");
				for (int f=nfaults; f<nfaults+Nfaults[nev]; f++) {
					for (int p=1; p<=(*eqk_aft)[f].np_di*(*eqk_aft)[f].np_st; p++) {
						for (int l=0; l<Nas; l++) {
							if ((*eqk_aft)[f].allslip_str) fprintf(fout,"%.5e\t",(*eqk_aft)[f].allslip_str[l][p]);
						}
						fprintf(fout,"\n");
					}
				}
			fclose(fout);
			nfaults+=Nfaults[nev];
		}
	}

	// Temporal evolution of afterslip.//
	nfaults=0;
	for (int nev=0; nev<NA; nev++){
		Teq=(*eqk_aft)[nfaults].t;
		t_afterslip=(*eqk_aft)[nfaults].ts;
		Nas=(*eqk_aft)[nfaults].nosnap;
		for (int i=0; i<Nas; i++) t_afterslip[i]-=Teq;	//since functions below start from t=0;
		for (int i=0; i<=*L; i++) (*times2)[i]-=Teq;	//since functions below start from t=0;

		splines_eqkfm(&eq_aft, Nas, Nfaults[nev], t_afterslip-1, (*times2)-1, *L+1, seed);
		//at this point eq_aft[f].allslip_str[l][p] contains cumulative slip on patch [p] at time times2[l].

		for (int f=0; f<Nfaults[nev]; f++) {

			for (int p=1; p<=eq_aft[f].np_di*eq_aft[f].np_st; p++) {
				for (int l=*L; l>=0; l--) {
					if ((*times2)[l]<0.0){
						//no afterslip before its start time:
						if (eq_aft[f].allslip_str) eq_aft[f].allslip_str[l][p]=0.0;
						if (eq_aft[f].allslip_dip) eq_aft[f].allslip_dip[l][p]=0.0;
						if (eq_aft[f].allslip_open) eq_aft[f].allslip_open[l][p]=0.0;
					}
					else{
						if (l>0 && (*times2)[l-1]>0.0){	//this is to avoid subtracting from element with t<Teq.
							if (eq_aft[f].allslip_str) eq_aft[f].allslip_str[l][p]-=eq_aft[f].allslip_str[l-1][p];
							if (eq_aft[f].allslip_dip) eq_aft[f].allslip_dip[l][p]-=eq_aft[f].allslip_dip[l-1][p];
							if (eq_aft[f].allslip_open) eq_aft[f].allslip_open[l][p]-=eq_aft[f].allslip_open[l-1][p];
						}
					}
				}

				//shift by one since eq_aft[f].allslip_str[l][p] must refer to time between times2[t+1] and times2[t].
				for (int l=0; l<*L; l++) {
					if (eq_aft[f].allslip_str) eq_aft[f].allslip_str[l][p]=eq_aft[f].allslip_str[l+1][p];
					if (eq_aft[f].allslip_dip) eq_aft[f].allslip_dip[l][p]=eq_aft[f].allslip_dip[l+1][p];
					if (eq_aft[f].allslip_open) eq_aft[f].allslip_open[l][p]=eq_aft[f].allslip_open[l+1][p];
				}
			 }
			eq_aft[f].tevol=NULL;
		}

		for (int i=0; i<Nas; i++) t_afterslip[i]+=Teq;	//revert shift from before;
		for (int i=0; i<=*L; i++) (*times2)[i]+=Teq;	//revert shift from before;
		nfaults+=Nfaults[nev];
		eq_aft+=Nfaults[nev];	//gets shifted every time so 0th element is the one of the new nev.
	}

	*eqk_aft=eq_aft-nfaults;	//shift it back.

	//printout TODO add str, dip, open. (also above).

	if (printout_splines){
		nfaults=0;
		for (int nev=0; nev<NA; nev++){
			sprintf(fname,"splines%d.dat",nev);
			fout=fopen(fname,"w");
			for (int l=0; l<=*L-1; l++) {
				fprintf(fout,"%.5e\t",(*times2)[l]);
			}
			fprintf(fout,"\n\t");
				for (int f=nfaults; f<nfaults+Nfaults[nev]; f++) {
					for (int p=1; p<=(*eqk_aft)[f].np_di*(*eqk_aft)[f].np_st; p++) {
						for (int l=0; l<=*L-1; l++) {
							if ((*eqk_aft)[f].allslip_str) fprintf(fout,"%.5e\t",(*eqk_aft)[f].allslip_str[l][p]);
						}
						fprintf(fout,"\n");
					}
				}
			fclose(fout);
			nfaults+=Nfaults[nev];
		}
	}

	return(err!=0);
}

int setup_afterslip_single_log(double t0, double t1, double *Cs, double *ts,
						 int Nfun, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2,
						 long *seed) {

//if splines==1, eqk_aft is substituted with more densily discretized version (L steps instead of Nas).
//Cs, ts, coefficients of temporal evolution functions (See Savage Parkfield paper). Nfun: no. of such funtions.
// eq_aft has indices: [0...Nas*Nfaults-1].
// logarithmic: flag indicating if stressing history is logarithmic. If not, use linear.

	int err=0;
	double TAU=200000, dtau=7000, timeTAU=183;	//todo allow to set from outside.
	double smallstepstime=12;
	double now, prev, norm, curr;
	double Tendaft;	//Time to which cumulative afterslip snapshot refers.
	struct eqkfm *eq_aft;

	double Teq;
	int nfaults;

	eq_aft= *eqk_aft;

	err=timesteps_omori(t0, t1, eqk_aft, NA, Nfaults, L, times2, smallstepstime, TAU, dtau, timeTAU);

	// Temporal evolution of afterslip.//
	nfaults=0;
	for (int nev=0; nev<NA; nev++){
		Tendaft=(*eqk_aft)[nfaults].ts[(*eqk_aft)[nfaults].nosnap-1];
		Teq=(*eqk_aft)[nfaults].t;

		//allocate tevol vector:
		(*eqk_aft)[nfaults].tevol=dvector(0,*L-1);

		norm=0.0;
		for (int i=0; i<Nfun; i++) norm+=Cs[i]*log(1+(Tendaft-Teq)/ts[i]);
		now=0.0;
		curr=0.0;

		for (int t=1; t<=*L; t++){
			if ((*times2)[t]<Teq){
				(*eqk_aft)[nfaults].tevol[t-1]= 0.0;
			}
			else{
				prev=curr;
				now= (*times2)[t]-Teq;
				curr=0.0;
				for (int i=0; i<Nfun; i++) curr+=Cs[i]*log(1+now/ts[i]);
				(*eqk_aft)[nfaults].tevol[t-1]= (curr-prev)/norm;
			}
		}

		//assign tevol vectors to all subfaults:
		for (int f=0; f<Nfaults[nev]; f++) {
			(*eqk_aft)[nfaults+f].tevol=(*eqk_aft)[nfaults].tevol;
		}

		nfaults+=Nfaults[nev];
	}

	return(err!=0);
}