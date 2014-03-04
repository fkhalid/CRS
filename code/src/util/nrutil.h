/******************************************************************************
  "nrutil.c - Numerical Recipes utility file"

  description:  numerical recipes for vector and matrix declaration

  (C) Copr. 1986-92 Numerical Recipes Software

******************************************************************************/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

/* Numerical Recipes standard error handler */
void nrerror(char error_text[]);

/* allocate a float vector with subscript range v[nl..nh] */
float *vector(long nl, long nh);

/* allocate an int vector with subscript range v[nl..nh] */
int *ivector(long nl, long nh);

/* allocate a short int vector with subscript range v[nl..nh] */
short int *sivector(long nl, long nh);

/* allocate an unsigned char vector with subscript range v[nl..nh] */
unsigned char *cvector(long nl, long nh);

/* allocate an unsigned long vector with subscript range v[nl..nh] */
unsigned long *lvector(long nl, long nh);

/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(long nl, long nh);

/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
float **matrix(long nrl, long nrh, long ncl, long nch);

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl, long nch);

/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
int **imatrix(long nrl, long nrh, long ncl, long nch);

/* allocate a short int matrix with subscript range m[nrl..nrh][ncl..nch] */
short int **simatrix(long nrl, long nrh, long ncl, long nch);

/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl);

/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);

/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

/* free a float vector allocated with vector() */
void free_vector(float *v, long nl, long nh);

/* free an int vector allocated with ivector() */
void free_ivector(int *v, long nl, long nh);

/* free an unsigned char vector allocated with cvector() */
void free_cvector(unsigned char *v, long nl, long nh);

/* free an unsigned long vector allocated with lvector() */
void free_lvector(unsigned long *v, long nl, long nh);

/* free a double vector allocated with dvector() */
void free_dvector(double *v, long nl, long nh);

/* free a float matrix allocated by matrix() */
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);

/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

/* free an int matrix allocated by imatrix() */
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

/* free a submatrix allocated by submatrix() */
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);

/* free a matrix allocated by convert_matrix() */
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);

/* free a float f3tensor allocated by f3tensor() */
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
