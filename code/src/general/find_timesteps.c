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

int findtimestepsomori(double te, double t0, double t1, double K, double p, double c, double *t, double *K0, int *L) {
/*
 * Calculates time steps with an increasing spacing, a function of the form dt~(t+c)^p is used: for p=1 and a logarithmic stressing history,
 * this gives equal stresses within each time step. (NB: t refers to the start of each aseismic event).
 *
 * Input:
 * 	te: earthquake time (or in general, start time of the aseismic event).
 * 	t0,t1: start and end time of the period for which time steps are calculated.
 *
 * 	K, p, c: parameters controlling the shape of the function used to estimate time steps: *
 *
 * Output:
 *  t is populated with the time steps (NB must be allocated beforehand).
 *  K0: value such that t_{i}=t_{i-1}+K(t+c-teq)^p. Ignored if NULL.
 *  L: number of time steps.
 */


	// [Fahad] Variables used for MPI

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double dfdt, dt, tnow;
	int N, j=0, err=0;

	dt=K*pow(t0+c-te,p);
	N=(int) (t1-t0)*(1.0/dt);			//N is always > than the needed number of elements since the first derivative is monotonically decreasing.

	tnow=t0;
	if (t) t[0]=t0;
	while (tnow<=t1){
		dt=K*pow(tnow+c-te,p);
		tnow+=dt;
		if (t) t[j+1]=tnow;
		j+=1;
	}

	if (j>1000) print_screen("\n ** Warning: findtimesteps.m produced %d time steps! **\n",j);

	*L=j;
	if (K0) *K0= K;	//so that t_{i}=t_{i-1}+K(t+c-teq)^p



	return err;
}
