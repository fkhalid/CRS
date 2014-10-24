/*
 * moreutils.h
 *
 *  Created on: Feb 9, 2012
 *      Author: camcat
 */

# include <math.h>

int closest_element(double *v, int N, double value, double toll);
int *nth_index(int i, int Ndim, int *dim);
void nrerrorsoft(char error_text[]);
void copy_matrix( double **m1, double ***m2, int a, int b);
void mysort(unsigned long n, double *old_arr, int **ind, double **arr);
char ***tmatrix(long nrl, long nrh, long ncl, long nch, long length);
void free_tmatrix(char ***m, long nrl, long nrh, long ncl, long nch, long length);
double **mtimesm3(double **m1, double **m2, double ***);
double * mtimesv(double **M, double *v, double *v2, int D1, int D2);
double norm(double *v1, int D);
double vdotv(double *v1, double *v2, int D);
void normv (double *v, int D);
void unitv (double *v, int D);
void statistics (double *, int, double *, double *);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,	long ndl, long ndh);
void intersect_lists(int *l1, int *l2, int **l3, int **, int **, int N1, int N2, int *N3);
void nearest_neighbours(int NP, int D1, int D2, int D3, int **nn);
void interp_nn(int NP, int D1, int D2, int D3, double *values, double **allvalues, int all6, int **);
double *** duplicate_d3tensor(double ***S, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double * duplicate_dvector(double *v, long nrl, long nrh);
double min_v(double *v, int N);
double max_v(double *v, int N);


