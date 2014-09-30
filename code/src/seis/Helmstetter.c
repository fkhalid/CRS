/*
 * Helmstetter.c
 *
 *  Created on: Nov 8, 2013
 *      Author: camcat
 *
 * Produces smoothed seismicity model as in Helmstetter (2007), using a gaussian kernel.
 * Catalog need to be previously declustered.
 */

#include "Helmstetter.h"

double *Helmstetter(double *xgrid, double *ygrid, double dx, double dy, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord){
/* weights: flags used for declustering (0/1 for excluded/selected events). if NULL, all events are selected with weight=1.
 * ord= 1,2: indicates if first o fsecond nnearest neighbour should be used.
 *
 */

	double d, w;
	double *dist=NULL;
	int *ind, no_ind;
	double *rate, *rate_tot;

	ind=ivector(1,Ngrid);
	rate=dvector(1,Ngrid);
	rate_tot=dvector(1,Ngrid);
	for (int i=1; i<=Ngrid; i++) rate_tot[i]=0.0;

	switch (ord){
		case 1:
			all_nearestneighbours(xs, ys, N, NULL, &dist);
			break;
		case 2:
			all_2ndnearestneighbours(xs, ys, N, NULL, &dist);
			break;
		default:
			print_screen("** Error:  illegal value for variable 'ord' in Helmstetter.c. \n**");
			print_logfile("** Error:  illegal value for variable 'ord' in Helmstetter.c. \n**");
			return NULL;
	}

	for (int eq=1; eq<=N; eq++){	//todo parallel.
		if (!weights || (weights[eq]>0.0)){
			d=fmax(dist[eq],err[eq]);
			find_gridpoints_exact(ygrid, xgrid, NULL, dx, dy, 0.0, Ngrid, Ngrid, ys[eq], xs[eq], d, 0.0, 0.0, 10000, &no_ind, ind, rate, 1, 0);
			w= (weights)? weights[eq] : 1.0;
			for (int i=1; i<=no_ind; i++) rate_tot[ind[i]]+=w*rate[i];
		}
	}

	return rate_tot;
}

double *Helmstetter_nonuni(double *xgrid, double *ygrid, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord){
/* weights: flags used for declustering (0/1 for excluded/selected events). if NULL, all events are selected with weight=1.
 * ord= 1,2: indicates if first o fsecond nnearest neighbour should be used.
 * does not use exact prop. inside grid point, but gaussian distr. based on distance from center point.
 *
 */

	double d, w;
	double *dist=NULL;
	int *ind, no_ind;
	double *rate, *rate_tot;

	ind=ivector(1,Ngrid);
	rate=dvector(1,Ngrid);
	rate_tot=dvector(1,Ngrid);
	for (int i=1; i<=Ngrid; i++) rate_tot[i]=0.0;

	switch (ord){
		case 1:
			all_nearestneighbours(xs, ys, N, NULL, &dist);
			break;
		case 2:
			all_2ndnearestneighbours(xs, ys, N, NULL, &dist);
			break;
		default:
			print_screen("** Error:  illegal value for variable 'ord' in Helmstetter.c. \n**");
			print_logfile("** Error:  illegal value for variable 'ord' in Helmstetter.c. \n**");
			return NULL;
	}

	for (int eq=1; eq<=N; eq++){	//todo parallel.
		if (!weights || (weights[eq]>0.0)){
			d=fmax(dist[eq],err[eq]);
			find_gridpoints(ygrid, xgrid, NULL, NULL, Ngrid, Ngrid, ys[eq], xs[eq], d, 0.0, 0.0001, 1, &no_ind, ind, rate, 1, 0);
			w= (weights)? weights[eq] : 1.0;
			for (int i=1; i<=no_ind; i++) rate_tot[ind[i]]+=w*rate[i];
		}
	}

	return rate_tot;
}

double *Helmstetter_cat(struct catalog cat, struct crust crst, double *weights, int ord) {

	double dx, dy;

	dx=(DEG2RAD*crst.dlon)*Re*cos(DEG2RAD*crst.lat0);
	dy=(DEG2RAD*crst.dlat)*Re;

	return Helmstetter(crst.x, crst.y, dx, dy, crst.nLon*crst.nLat, cat.x0, cat.y0, cat.err, weights, cat.Z, ord);

}

double *fit_depth(double *zgrid, double dz, int Ngrid, double *zs, double *err, double *weights, int N){
/* uses values in err as st.dev. of depth.
 * weights: array contain weights of earthquakes; if NULL, all weights will be set to 1.
 */


	double *prob0, *prob;
	double probCum;
	double rz, w;

	prob0=dvector(1,Ngrid);
	prob=dvector(1,Ngrid);

	for (int n=1; n<=Ngrid; n++) prob[n]=0;


	for (int eq=1; eq<=N; eq++){	//todo parallel.
		probCum=0;
		for (int n=1; n<=Ngrid; n++){
			rz=fabs(zgrid[n]-zs[eq]);

			prob0[n]=exact_prob_1d(rz, dz, err[eq]);
		}
		w=(weights)? weights[eq] : 1.0;
		for (int n=1; n<=Ngrid; n++) prob[n]+=w*prob0[n];
	}

	free_dvector(prob0,1,Ngrid);

	return prob;
}


