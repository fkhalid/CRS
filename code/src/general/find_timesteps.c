/*
 * findtimesteps.c
 *
 *  Created on: Mar 27, 2012
 *      Author: camcat
 */


#include "find_timesteps.h"

int findtimestepsomori(double te, double t0,double t1, double tstart, double tend, double tau0, double Dtau, double p, double c, double *t, double *K_over_tau0, int *L){
//te is earthquake time. tstart, tend are the start and end time of the measured value tau0. t0,t1 are the start and end time of the period for which time steps are calculated.
	double dfdt, dt;
	int N, j=0, err=0;

	*K_over_tau0=1/(pow(c+tend-te,1.0-p)-pow(c+tstart-te,1.0-p));

	dfdt=(1.0-p)*(*K_over_tau0)*tau0*pow(t0+c-te,-p);
	dt=Dtau*(1.0/dfdt);
	N=(int) (t1-t0)*(1.0/dt);			//N is always > than the needed number of elements since the first derivative is monotonically decreasing.

	t[0]=t0;
	while (t[j]<=t1){
		dfdt=(1.0-p)*(*K_over_tau0)*tau0*pow(t[j]+c-te,-p);
		dt=Dtau*(1.0/dfdt);
		t[j+1]=t[j]+dt;
		j+=1;
		if (j==(*L)) {
			printf("** Error: L too small in findtimestepsomori.c. Exiting. **\n");
			if (flog){
				fprintf(flog,"** Error: L too small in findtimestepsomori.c. Exiting. **\n");
				fflush(flog);
			}
			err=1;
			break;
		}
	}

	if (j>1000) printf("\n ** Warning: findtimesteps.m produced %d time steps! **\n",j);
	*L=j;
	return err;

}

void findtimestepslog(double t0,double t1, double tau0, double Dtau, double w, double **times, int *L){

	double dfdt, dt;
	double *t;
	int N, j=0;

//	*K_over_tau=1/(pow(c+tend,1.0-p)-pow(c+tstart,1.0-p));

	dfdt=tau0*w/(w*t0+1);
	dt=Dtau*(1.0/dfdt);
	N=(int) (t1-t0)*(1.0/dt);			//N is always > than the needed number of elements since the first derivative is monotonically decreasing.

	t=dvector(0,N-1);
	t[0]=t0;

	while (t[j]<=t1){
		dfdt=tau0*w/(w*t[j]+1);
		dt=Dtau*(1.0/dfdt);
		t[j+1]=t[j]+dt;
		j+=1;
	}

	if (j>10000) printf("\n ** Warning: findtimesteps.m produced more than 10000 time steps! **\n");

	*times=t;
	*L=j-1;
}

void tevolomori(double te, double *times, double Kotau, double p, double c, double *tevol, int L){
//te = earthquake time

	double dum1, dum2;
	FILE *fout;

	fout =fopen("/home/des/camcat/temp/times.dat","w");
	if (times[0]+c-te==0) printf("\n***Achtung!*** times[0]+c-te==0\n\n");
	dum2=pow(times[0]+c-te,1.0-p);
	for (int y=0; y<=L; y++){
		dum1=dum2;
		dum2=pow(times[y+1]+c-te,1.0-p);
		tevol[y]=Kotau*(dum2-dum1);
		fprintf(fout, "%lf\t%lf\n",times[y],tevol[y]);
	}

	fclose(fout);
}
