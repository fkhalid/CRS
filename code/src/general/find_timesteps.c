/*
 * findtimesteps.c
 *
 *  Created on: Mar 27, 2012
 *      Author: camcat
 */

#include "find_timesteps.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int findtimestepsomori(double te, double t0, double t1, double tstart, double tend,
					   double tau0, double Dtau, double p, double c, double *t,
					   double *K, int *L) {
//te is earthquake time. tstart, tend are the start and end time of the measured value tau0. t0,t1 are the start and end time of the period for which time steps are calculated.
	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double dfdt, dt;
	int N, j=0, err=0;
	double K_over_tau;

	K_over_tau=1/(pow(c+tend-te,1.0-p)-pow(c+tstart-te,1.0-p));

	dfdt=(1.0-p)*K_over_tau*tau0*pow(t0+c-te,-p);
	dt=Dtau*(1.0/dfdt);
	N=(int) (t1-t0)*(1.0/dt);			//N is always > than the needed number of elements since the first derivative is monotonically decreasing.

	t[0]=t0;
	while (t[j]<=t1){
		dfdt=(1.0-p)*K_over_tau*tau0*pow(t[j]+c-te,-p);
		dt=Dtau*(1.0/dfdt);
		t[j+1]=t[j]+dt;
		j+=1;
		if (j==(*L)) {
			print_screen("** Error: L too small in findtimestepsomori.c. Exiting. **\n");
			print_logfile("** Error: L too small in findtimestepsomori.c. Exiting. **\n");
			err=1;
			break;
		}
	}

	if (j>1000) print_screen("\n ** Warning: findtimesteps.m produced %d time steps! **\n",j);

	*L=j;
	if (K) *K= Dtau/((1.0-p)*K_over_tau*tau0);	//so that t_{i}=t_{i-1}+K(t+c-teq)^p

	return err;
}
