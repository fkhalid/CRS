/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.c.  Do not confuse this file with the same-named
   file nrutil.c that is supplied in the same subdirectory or archive
   as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil_newnames.h"

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

//float *vector(long nl, long nh)
///* allocate a float vector with subscript range v[nl..nh] */
//{
//	float *v;
//
//	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
//	if (!v) nrerror("allocation failure in vector()");
//	return v-nl+NR_END;
//}

int *iarray(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in iarray()");
	return v-nl+NR_END;
}

//unsigned char *cvector(long nl, long nh)
///* allocate an unsigned char vector with subscript range v[nl..nh] */
//{
//	unsigned char *v;
//
//	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
//	if (!v) nrerror("allocation failure in cvector()");
//	return v-nl+NR_END;
//}
//
unsigned long *larray(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in larray()");
	return v-nl+NR_END;
}

double *darray(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in darray()");
	return v-nl+NR_END;
}

float **f2array(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **d2array(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **i2array(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

// float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
// 	long newrl, long newcl)
// /* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
// {
// 	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
// 	float **m;
// 
// 	/* allocate array of pointers to rows */
// 	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
// 	if (!m) nrerror("allocation failure in submatrix()");
// 	m += NR_END;
// 	m -= newrl;
// 
// 	/* set pointers to rows */
// 	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
// 
// 	/* return pointer to array of pointers to rows */
// 	return m;
// }
// 
// float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
// /* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
// declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
// and ncol=nch-ncl+1. The routine should be called with the address
// &a[0][0] as the first argument. */
// {
// 	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
// 	float **m;
// 
// 	/* allocate pointers to rows */
// 	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
// 	if (!m) nrerror("allocation failure in convert_matrix()");
// 	m += NR_END;
// 	m -= nrl;
// 
// 	/* set pointers to rows */
// 	m[nrl]=a-ncl;
// 	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
// 	/* return pointer to array of pointers to rows */
// 	return m;
// }
// 
float ***f3array(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3array()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3array()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3array()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_iarray(int *v, long nl, long nh)
/* free an int vector allocated with iarray() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

//void free_cvector(unsigned char *v, long nl, long nh)
///* free an unsigned char vector allocated with cvector() */
//{
//	free((FREE_ARG) (v+nl-NR_END));
//}

void free_larray(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with larray() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_darray(double *v, long nl, long nh)
/* free a double vector allocated with darray() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_f2array(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by f2array() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_d2array(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by d2array() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_i2array(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by i2array() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

//void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
///* free a submatrix allocated by submatrix() */
//{
//	free((FREE_ARG) (b+nrl-NR_END));
//}
//
//void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
///* free a matrix allocated by convert_matrix() */
//{
//	free((FREE_ARG) (b+nrl-NR_END));
//}
//
void free_f3array(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3array allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

float ran1(long *idum)
{
/*	"Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added safeguards.
	Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).
	Call with idum a negative integer to initialize: thereafter, do not alter idum between successive deviates in a sequence.
	RNMX should approximate the largest floating value that is less than 1.
*/

	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {			//Initialize.
		if (-(*idum) < 1) *idum=1;		//Be sure to prevent idum = 0.
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {		//Load the shuffle table (after 8 warm-ups).
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;						//Start here when not initializing.
	*idum=IA*(*idum-k*IQ)-IR*k;			//Compute idum=(IA*idum) % IM without overflows by Schrage's method.
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;							//Will be in the range 0..NTAB-1.
	iy=iv[j];							//Output previously stored value and refill the shuffle table.
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;	//Because users don't expect endpoint values.
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software )!0,". */
