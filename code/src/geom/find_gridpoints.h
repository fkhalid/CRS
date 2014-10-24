/*
 * find_gridpoints.h
 *
 *  Created on: Dec 22, 2011
 *      Author: camcat
 */


#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../util/moreutil.h"
#include "../util/nrutil.h"

#define PI (3.141592653589793)

int find_gridpoints(double *ys, double *xs, double *dAs, double *depths, int N, double y, double x, double SD, double Depth, double SDd, int cut_sd, int *ngridj, int **ngridpointj, double **weightsj, int inside, int d3);
int find_gridpoints_d(double *ys, double *xs, double *depths, int *already_selected, int Nsel0, int N, double y_eq, double x_eq, double Depth, double m, double dDCFS, int *ngridj, int **ngridpointj);
int find_gridpoints_exact(double *ys, double *xs, double *depths, double dx, double dy, double dz, int N, int Nselmax, double y, double x, double SD, double Depth, double SDd,
		int cut_sd, int *ngridj, int *ngridpointj, double *weightsj, int inside, int d3);
double exact_prob_1d(double r, double dr, double sd);
double exact_prob(double rx, double ry, double rz, double dx, double dy, double dz, double sdx, double sdy, double sdz, int d3);
int all_nearestneighbours(double *x, double *y, int N, int **pts, double **dist);
int all_2ndnearestneighbours(double *x, double *y, int N, int **pts, double **dist);
