/*
 * find_gridpoints.c
 *
 *  Created on: Dec 22, 2011
 *      Author: camcat
 */

#include "find_gridpoints.h"


int find_gridpoints(double *ys, double *xs, double *dAs, double *depths, int N, int Nselmax, double y, double x, double SD, double Depth, double SDd,
		int cut_sd, int *ngridj0, int *ngridpointj, double *weightsj, int inside, int d3){

	double r, rmin=1e30, dz, probCum, *prob;
	int p,p2;
	int K, Kd;
	int closestp;
	int ngridj;
	double y1, y2, x1, x2, D1, D2, A, Vfrac=1;

	prob=dvector(0,N+1);
	if (!dAs) inside=1; 	//can't calculate total area, so assume it's all inside.


	//todo delete
	fprintf(flog, "find_gridpoints:\n \t ys=%.3lf, %.3lf, %.3lf, ...\n", ys[1], ys[2], ys[N]);
	fprintf(flog, "\t xs=%.3lf, %.3lf, %.3lf, ...\n", xs[1], xs[2], xs[N]);
	fprintf(flog, "\t dAs=%.3lf, %.3lf, %.3lf, ...\n", dAs[1], dAs[2], dAs[N]);
	fprintf(flog, "\t zs=%.3lf, %.3lf, %.3lf, ...\n", depths[1], depths[2], depths[N]);
	fprintf(flog, "\t N=%d, Nselmax=%d\n", N, Nselmax);
	fprintf(flog, "\t y=%.3lf, x=%.3lf, z=%.3lf\n", x, y, Depth);
	fprintf(flog, "\t SD=%.3lf, SDd=%.3lf\n", SD, SDd);
	fprintf(flog, "\t cut_sd=%d, inside=%d, d3=%d\n", cut_sd, inside, d3);

	ngridj=0;
	probCum=0;

// K, Kd determine cutoff radius/depth.
	K=cut_sd;
	Kd=cut_sd;

	y1=y-K*SD;
	y2=y+K*SD;
	x1=x-K*SD;
	x2=x+K*SD;
	D1=Depth-Kd*SDd;
	D2=Depth+Kd*SDd;

	A=0;

	for (p=1; p<=N; p++){
		r= (d3)? sqrt(pow(ys[p]-y,2)+pow(xs[p]-x,2)+pow(depths[p]-Depth,2)) : sqrt(pow(ys[p]-y,2)+pow(xs[p]-x,2));
		if (r<rmin){
			rmin=r;
			closestp=p;
		}
		if (ys[p]>=y1 && ys[p]<=y2 && xs[p]>=x1 && xs[p]<=x2 && (!d3 || (depths[p]>=D1 && depths[p]<=D2))){
			dz= (d3) ? depths[p]-Depth : 0.0;

			if (r<=K*SD && dz<=Kd*SDd)
			{
				ngridj+=1;
				if (ngridj>Nselmax){
					if (verbose_level) printf("*Error: *ngridj>Nselmax in find_gridpoints.c - need to choose larger value for Nselmax. Exiting. **\n");
					if (flog){
						fprintf(flog, "*Error: *ngridj>Nselmax in find_gridpoints.c - need to choose larger value for Nselmax. Exiting. **\n");
						fflush(flog);

					}
					return(1);
				}
				ngridpointj[ngridj]=p;
				prob[ngridj]= (d3)? exp(-pow(r,2)/(2*pow(SD,2)))*exp(-pow(dz,2)/(2*pow(SDd,2))) : exp(-pow(r,2)/(2*pow(SD,2)));
				probCum+=prob[ngridj];
				if (dAs) A+= dAs[p];
			}
		}
	}

// K*SD is cutoff radius (gaussian would imply inf points), so that total area considered is pi*(K*SD). Vfrac is the fraction of area inside grid (if event is located outside grid).
// Vfrac is fraction of area inside selected area.
	switch (inside) {
	case 0:
	  Vfrac=fabs(A/(pi*(K*SD)*(K*SD)));
	  break;
	case 1:
	  Vfrac=1;
	  break;
	default:
	  if (verbose_level>1) printf("Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
		if (flog) {
			fprintf(flog, "Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
			fflush(flog);
		}
	  return 1;
	  break;
	}

	if (Vfrac>1.1) {
		if (verbose_level>1) printf("Warning: Vfrac>1 (%lf)\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
		if (flog) {
			fprintf(flog,"Warning: Vfrac>1 (%lf) \n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
			fflush(flog);
		}
	}
	if (Vfrac>1) Vfrac=1;

	for (p2=1; p2<=ngridj; p2++) weightsj[p2]=Vfrac*prob[p2]/probCum;
	weightsj[0]=1-Vfrac;

	//if no point is selected, select nearest point:
	if (ngridj==0){
		ngridj=1;
		ngridpointj[1]=closestp;
		weightsj[1]=1;
	}

	if (ngridj0) *ngridj0=ngridj;
	free_dvector(prob, 0, N+1);
	return(0);
}

int find_gridpoints_d(double *ys, double *xs, double *depths, int *already_selected, int Nsel0, int N, double y_eq, double x_eq, double Depth, double m,
		double dDCFS, int *ngridj, int **ngridpointj){
// xs, ys: eastwards, northwards.
// this function finds gridpoints within given distance(without weighting them as in find_gridpoints).
// already_selected are the points selected as possible location of earthquake (should all be included). There are Nsel0 of them.

	int *points_temp;
	double r, rmin=1e30;
	double y1, y2, x1, x2, D1, D2;
	double k=0.85, R;	//k found empirically looking at distrib. of DCFS(r,m). dcfs_min is threshold for coulomb stress value, in Pa (used to select points: points that are furhter than distance corresp. to this are ignored).
	int counter=1, closestp;

	points_temp=ivector(1,N);
	*ngridj=0;

	R=pow(k*pow(10.0,3.0*m/2.0)/dDCFS,1.0/3.0);
	y1=y_eq-R;
	y2=y_eq+R;
	x1=x_eq-R;
	x2=x_eq+R;
	D1=fmax(Depth-R, 0.0);
	D2=Depth+R;

	for (int p=1; p<=N; p++){
		if (counter<=Nsel0 && already_selected[counter]==p){
			counter+=1;
			*ngridj+=1;
			points_temp[*ngridj]=p;
		}
		else {
			r=sqrt(pow(ys[p]-y_eq,2)+pow(xs[p]-x_eq,2)+pow(depths[p]-Depth,2));
			if (r<rmin){
				rmin=r;
				closestp=p;
			}
			if (ys[p]>=y1 && ys[p]<=y2 && xs[p]>=x1 && xs[p]<=x2 && depths[p]>=D1 && depths[p]<=D2){
				if (r<=R){
					*ngridj+=1;
					points_temp[*ngridj]=p;
				}
			}
		}
	}

	if (*ngridj==0){
		*ngridj=1;
		*ngridpointj=ivector(1,(*ngridj));
		(*ngridpointj)[1]=closestp;
	}

	else{
		*ngridpointj=ivector(1,(*ngridj));
		for (int p2=1; p2<=*ngridj; p2++) (*ngridpointj)[p2]= points_temp[p2];
	}
	free_ivector(points_temp,1,N);

	return(0);
}

int find_gridpoints_exact(double *ys, double *xs, double *depths, double dx, double dy, double dz, int N, int Nselmax, double y, double x,
		double SD, double Depth, double SDd, int cut_sd, int *ngridj, int *ngridpointj, double *weightsj, int inside, int d3){

/* Instead of simply using center point, integrates over each cell.
 * d3= use 3d distance (as opposed to horizontal).
 * if inside==1, sum of weights is 1; otherwise, is can be smaller than 1 if part of the gaussian is outside of the domain of xs, ys, depths.
 *
 * ngridj=no. of points selected. Ignored if NULL is passed.
 */

	double r, rmin=1e30, probCum, prob[N+1];
	double rx, ry, rz;
	int p,p2;
	int K, Kd;
	int closestp, ngridj_int=0;
	double y1, y2, x1, x2, D1, D2, A, Vfrac=1;

	if (ngridj) *ngridj=0;
	probCum=0;

// K, Kd determine cutoff radius/depth.
	K=cut_sd;
	Kd=cut_sd;

	y1=y-K*SD;
	y2=y+K*SD;
	x1=x-K*SD;
	x2=x+K*SD;
	D1=Depth-Kd*SDd;
	D2=Depth+Kd*SDd;

	for (p=1; p<=N; p++){
		rx=xs[p]-x;
		ry=ys[p]-y;
		rz= (d3)? depths[p]-Depth : 0.0;
		r= sqrt(pow(ry,2)+pow(rx,2)+pow(rz,2));
		if (r<rmin){
			rmin=r;
			closestp=p;
		}
		if (ys[p]>=y1 && ys[p]<=y2 && xs[p]>=x1 && xs[p]<=x2 && (!d3 || (depths[p]>=D1 && depths[p]<=D2))){
			if (r<=K*SD && (!d3 || rz<=Kd*SDd)){
				ngridj_int+=1;
				if (ngridj_int>Nselmax){
					if (verbose_level) printf("*Error: *ngridj>Nselmax in find_gridpoints.c - need to choose larger value for Nselmax. Exiting. **\n");
					if (flog){
						fprintf(flog, "*Error: *ngridj>Nselmax in find_gridpoints.c - need to choose larger value for Nselmax. Exiting. **\n");
						fflush(flog);

					}
					return(1);
				}
				ngridpointj[ngridj_int]=p;
				prob[ngridj_int]= exact_prob(rx,ry,rz,dx, dy,dz,SD, SD, SDd, d3);
				probCum+=prob[ngridj_int];
			}
		}
	}

	A=ngridj_int*dx*dy;

// K*SD is cutoff radius (gaussian would imply inf points), so that total area considered is pi*(K*SD). Vfrac is the fraction of area inside grid (if event is located outside grid).
// Vfrac is fraction of area inside selected area.

	switch (inside) {
		case 0:
		  Vfrac=fabs(A/(pi*(K*SD)*(K*SD)));
		  break;
		case 1:
		  Vfrac=1;
		  break;
		default:
	  if (verbose_level>1) printf("Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
		if (flog) {
			fprintf(flog, "Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
			fflush(flog);
		}
	  return 1;
		  break;
	}

	if (Vfrac>1.1) {
		if (verbose_level>1) printf("Warning: Vfrac>1 (%lf)\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
		if (flog) {
			fprintf(flog,"Warning: Vfrac>1 (%lf) \n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
			fflush(flog);
		}
	}
	if (Vfrac>1) Vfrac=1;

	for (p2=1; p2<=ngridj_int; p2++) weightsj[p2]=Vfrac*prob[p2]/probCum;
	weightsj[0]=1-Vfrac;

	if (ngridj_int==0){
		ngridj_int=1;
		ngridpointj[1]=closestp;
		weightsj[1]=1;
	}

	if (ngridj) *ngridj=ngridj_int;

	return(0);
}

double exact_prob_1d(double r, double dr, double sd){
/* r= distance of cell center to central point of gaussian; \
 * dr= cell size;
 * sd= st.dev. of gaussian.
 */

	double 	rmin = r-0.5*dr,	\
			rmax = r+0.5*dr;

	return (erf(rmax/(sqrt(2.0)*sd))-erf(rmin/(sqrt(2.0)*sd)))/2.0;

}


double exact_prob(double rx, double ry, double rz, double dx, double dy, double dz, double sdx, double sdy, double sdz, int d3){
	/* fod each dimension:
	 * rX= distance of cell center to central point of gaussian; \
	 * dX= cell size;
	 * sdX= st.dev. of gaussian.
	 * d3= flag to indicate is 3D distance should be used (otherwise, ignore vertical coordinate).
	 */

	double 	xmin = rx-0.5*dx,	\
			xmax = rx+0.5*dx,	\
			ymin = ry-0.5*dy,	\
			ymax = ry+0.5*dy,	\
			zmin = rz-0.5*dz,	\
			zmax = rz+0.5*dz;

	double Ix, Iy, Iz, I;

	Ix=erf(xmax/(sqrt(2.0)*sdx))-erf(xmin/(sqrt(2.0)*sdx));
	Iy=erf(ymax/(sqrt(2.0)*sdy))-erf(ymin/(sqrt(2.0)*sdy));
	if (d3) Iz=erf(zmax/(sqrt(2.0)*sdz))-erf(zmin/(sqrt(2.0)*sdz));

	I= (d3)? Ix*Iy*Iz/(8*dx*dy*dz) : Ix*Iy/(4*dx*dy);

	return I;

}

int all_nearestneighbours(double *x, double *y, int N, int **pts, double **dist){
	/* x, y: coordinates (indices: [1...N]);
	 * pts: index of nearest neighbour;
	 * dist: distance to nearest neighbour;
	 * pts, dist are pointers to 1D arrays. If NULL, ignored; if they point to NULL, memory will be allocated. Otherwise, arrays of the correct size should be passed.
	 */

	double *xs=NULL;	//sorted copies of x,y, (will be sorted).
	int *ind=NULL;
	int *x_ind=NULL;
	int *x_order;
	double d, *dmin= NULL;
	double y_indp;
	int x_toofar, x_toofarahead, x_toofarbehind;
	int indp0, indp;
	int n_op=0;	//no. of times distance is computed.

	if (pts) ind=*pts;
	if (dist) dmin=*dist;
	if (!dmin) {
		dmin=dvector(1,N);
		if (dist) *dist=dmin;
	}
	if (!ind) {
		ind=ivector(1,N);
		if (pts) *pts=ind;
	}
	for (int i=1; i<=N; i++) dmin[i]=1e30;

	// sort element by x;

	mysort(N, x, &x_ind, &xs);
	x_order=ivector(1,N);
	for (int i=1; i<=N; i++) x_order[x_ind[i]]=i;

	for (int i=1; i<=N; i++){
		indp0=indp=x_order[i];
		x_toofar=x_toofarahead=x_toofarbehind=0;
		//search in x direction.
		while (!x_toofar){
			indp=(indp>indp0)? 2*indp0-indp : 2*indp0-indp+1;
			if (indp>indp0){
				if (x_toofarahead) continue;
				if (indp>N || fabs(xs[indp]-x[i])>dmin[i]) {
					x_toofarahead=1;
					if (x_toofarbehind) x_toofar=1;
					continue;
				}
			}
			else{
				if (x_toofarbehind) continue;
				if (indp<=0 || fabs(xs[indp]-x[i])>dmin[i]) {
					x_toofarbehind=1;
					if (x_toofarahead) x_toofar=1;
					continue;
				}
			}
			y_indp= y[x_ind[indp]];
			if (fabs(y_indp-y[i])> dmin[i]) continue;
			d=sqrt(pow(xs[indp]-x[i],2)+pow(y_indp-y[i],2));
			n_op+=1;
			if (d<dmin[i]){
				dmin[i]=d;
				ind[i]=x_ind[indp];
			}
			if (d<dmin[x_ind[indp]]){
				dmin[x_ind[indp]]=d;
				ind[x_ind[indp]]=i;
			}
		}
	}

	return n_op;
}


int all_2ndnearestneighbours(double *x, double *y, int N, int **pts, double **dist){
	/* x, y: coordinates (indices: [1...N]);
	 * pts: index of nearest neighbour;
	 * dist: distance to nearest neighbour;
	 * pts, dist are pointers to 1D arrays. If NULL, ignored; if they point to NULL, memory will be allocated. Otherwise, arrays of the correct size should be passed.
	 */

	double *xs=NULL;	//sorted copies of x,y, (will be sorted).
	int **ind=NULL;
	int *x_ind=NULL;
	int *x_order;
	double d, **dmin;
	double y_indp;
	int x_toofar, x_toofarahead, x_toofarbehind;
	int indp0, indp;
	int n_op=0;	//no. of times distance is computed.

	dmin=dmatrix(1,2,1,N);
	ind=imatrix(1,2,1,N);

	for (int i=1; i<=N; i++) dmin[1][i]=dmin[2][i]=1e30;

	// sort element by x;
	mysort(N, x, &x_ind, &xs);
	x_order=ivector(1,N);
	for (int i=1; i<=N; i++) x_order[x_ind[i]]=i;

	for (int i=1; i<=N; i++){
		indp0=indp=x_order[i];
		x_toofar=x_toofarahead=x_toofarbehind=0;
		//search in x direction.
		while (!x_toofar){
			indp=(indp>indp0)? 2*indp0-indp : 2*indp0-indp+1;
			if (indp>indp0){
				if (x_toofarahead) continue;
				if (indp>N || fabs(xs[indp]-x[i])>dmin[2][i]) {
					x_toofarahead=1;
					if (x_toofarbehind) x_toofar=1;
					continue;
				}
			}
			else{
				if (x_toofarbehind) continue;
				if (indp<=0 || fabs(xs[indp]-x[i])>dmin[2][i]) {
					x_toofarbehind=1;
					if (x_toofarahead) x_toofar=1;
					continue;
				}
			}
			y_indp= y[x_ind[indp]];
			if (fabs(y_indp-y[i])> dmin[2][i]) continue;
			d=sqrt(pow(xs[indp]-x[i],2)+pow(y_indp-y[i],2));
			n_op+=1;
			if (d<dmin[2][i] && x_ind[indp]!=ind[1][i]){
				if (d<dmin[1][i]){
					dmin[2][i]=dmin[1][i];
					ind[2][i]=ind[1][i];
					dmin[1][i]=d;
					ind[1][i]=x_ind[indp];
				}
				else {
					dmin[2][i]=d;
					ind[2][i]=x_ind[indp];
				}
			}

			if (d<dmin[2][x_ind[indp]] && ind[1][x_ind[indp]]!=i){
				if (d<dmin[1][x_ind[indp]]){
					dmin[2][x_ind[indp]]=dmin[1][x_ind[indp]];
					ind[2][x_ind[indp]]=ind[1][x_ind[indp]];
					dmin[1][x_ind[indp]]=d;
					ind[1][x_ind[indp]]=i;
				}
				else {
					dmin[2][x_ind[indp]]=d;
					ind[2][x_ind[indp]]=i;
				}
			}
		}
	}

	if (pts) *pts=ind[2];
	if (dist) *dist=dmin[2];

	return n_op;
}
