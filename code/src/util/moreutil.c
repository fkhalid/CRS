
/*   Copyright (C) 2015 by Camilla Cattania and Fahad Khalid.
 *
 *   This file is part of CRS.
 *
 *   CRS is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CRS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CRS.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "moreutil.h"

#include <math.h>
#include <nrutil.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define NSTACK 50
#define FREE_ARG char*
#define NR_END 1

int closest_element(double *v, int N, double value, double toll){
	/*
	 * Returns index of the element in v closes to value; returns -1 if the difference is larger than toll.
	 * v range [0...N-1]
	 */

	int i;
	double mindist=1e30;

	for (int n=0; n<N; n++){
		if(fabs(v[n]-value)<mindist){
			mindist=fabs(v[n]-value);
			i=n;
		}
	}

	if (mindist>toll) return -1;
	else return i;
}

int *nth_index(int i, int Ndim, int *dim){
//given linear index and N dimensions with indices [1...dim[0]], [1...dim[1]], ..., [1...dim[Ndim-1]], returns array with indices of each dimension.

	int f=1;
	int *res=ivector(0,Ndim-1);

	for (int d=0; d<Ndim; d++) {
		res[d]=((i-1)%(f*dim[d]))/f+1;
		f*=dim[d];
	}

	return res;
}

// error handling w/o exiting code (hence "soft")
void nrerrorsoft(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...");
	fprintf(stderr,"%s\n",error_text);
}

void copy_matrix( double **m1, double ***m2, int a, int b){
	/* m1: original matrix; *m2: new matrix. indices: [1...a, 1...b]
	 * if *m2==NULL, memory allocated. Otherwise, m2 must have correct no. of elements.
	 */

	if (!(*m2)) *m2=dmatrix(1,a,1,b);

	for (int ns=1; ns<=a; ns++){
		for (int n=1; n<=b; n++) (*m2)[ns][n]=m1[ns][n];
	}

	return;
}

void copy_vector(double *m1, double **m2, int a){
	/* m1: original vector; *m2: new vector. indices: [1...a]
	 * if *m2==NULL, memory allocated. Otherwise, m2 must have correct no. of elements.
	 */

	if (!(*m2)) *m2=dvector(1,a);

	for (int ns=1; ns<=a; ns++){
		(*m2)[ns]=m1[ns];
	}

	return;
}

void mysort(unsigned long n, double *old_arr, int **ind_out, double **arr_out){
	/*
	 * Like numerical recipes, but returns array containing old element index, and uses doubles.
	 * ind_out, arr_out are pointers to 1d arrays. If they point to NULL, memory will be allocated. Otherwise, arrays of the correct size should be passed.
	 */

	unsigned long i,ir=n,j,k,l=1,*istack;
	int jstack=0, M=7, ii;
	double a,temp;
	int *ind= (ind_out)? *ind_out : NULL;
	double *arr= *arr_out;

	if (!ind) ind=ivector(1,n);
	if (!arr) arr=dvector(1,n);

	for (i=1; i<=n; i++){
		ind[i]=i;
		arr[i]=old_arr[i];
	}

	istack=lvector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				ii=ind[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					ind[i+1]=ind[i];
				}
				arr[i+1]=a;
				ind[i+1]=ii;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			SWAP(ind[k],ind[l+1])
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
				SWAP(ind[l],ind[ir])
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
				SWAP(ind[l+1],ind[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
				SWAP(ind[l],ind[l+1])
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			ii=ind[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
				SWAP(ind[i],ind[j]);
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			ind[l+1]=ind[j];
			ind[j]=ii;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_lvector(istack,1,NSTACK);

	if (ind_out) *ind_out=ind;
	if (arr_out) *arr_out=arr;

	return;
}

int **imatrix_firstlevel(long nrh)
/* allocate a int matrix with subscript range m[0..nrh]; columns will be allocated later. */
{
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrh+1)*sizeof(int*)));

	/* set rows to NULL (so that later it will be recognized that memory should be allocated).*/
	for (int i=0; i<=nrh; i++) m[i]=NULL;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix_firstlevel(long nrh)
/* allocate a double matrix with subscript range m[nrl..nrh];  columns will be allocated later. */
{
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrh+1)*sizeof(double*)));

	/* set rows to NULL (so that later it will be recognized that memory should be allocated).*/
	for (int i=0; i<=nrh; i++) m[i]=NULL;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_imatrix_firstlevel(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix_firstlevel() */
{
	for (int i=nrl; i<=nrh; i++) if (m[i]) free_ivector(m[i], ncl, nch);
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix_firstlevel(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by dmatrix_firstlevel() */
{
	for (int i=nrl; i<=nrh; i++) if (m[i]) free_dvector(m[i], ncl, nch);
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_tmatrix(char ***m, long nrl, long nrh, long ncl, long nch, long length)
/* free a char* matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl][ncl]+length-1-NR_END));
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

double **mtimesm3(double **m1, double **m2, double ***m30){
//if m30==NULL, m2 will be ignored and new memory wll be allocated; otherwise, m30 is filled (NB memory needs to have been allocated previously!).

//indices: [1...3].

	double **m3;

	m3= (m30) ? *m30 : dmatrix(1,3,1,3);

	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++){
			m3[i][j]=0;
			for (int k=1; k<=3; k++){
				m3[i][j]+=m1[i][k]*m2[k][j];
			}
		}
	}

	return m3;
}

double *mtimesv(double **M, double *v, double *v2, int D1, int D2){
//if v2==NULL, v2 will be ignored and new memory wll be allocated; otherwise, vector v2 is filled (NB memory needs to have been allocated previously!).

	int i, j;
	double temp;
	double *v2int;

	v2int=(v2)? v2 : dvector(1,D2);

	for (i=1; i<=D2; i++){
		temp=0;
		for (j=1;j<=D1; j++){
			temp+=M[i][j]*v[j];
		}
	v2int[i]=temp;
	}
	return v2int;
}

double norm(double *v1, int D){
//calculates vector norm.

	double r2;
	r2=vdotv(v1,v1,D);

	return sqrt(r2);
}

double vdotv(double *v1, double *v2, int D){
//dot product.

	int k;
	double temp2=0;
	for (k=1; k<=D; k++){
		temp2+=v1[k]*v2[k];
	}
	return temp2;
}

void normv (double *v, int D){
// normalizes vector v (so that the sum of its elements is 1).

	double sum=0.0;
	for (int k=1; k<=D; k++) sum+=v[k];
	for (int k=1; k<=D; k++) v[k]*=(1.0/sum);

	return;
}

void unitv (double *v, int D){
// normalizes vector v (so that its length is 1).

	double length2;
	length2=vdotv(v, v, D);
	for (int k=1; k<=D; k++){
			v[k]=v[k]/pow(length2,0.5);
		}
	return;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
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

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void nearest_neighbours(int NP, int D1, int D2, int D3, int **nn){
// points change along D1, then along D2, then along D3 (DX is no of points in each dimension).
	int P[4], P2[4];
	int Ds[4]={0, D1, D2, D3};
	int D1D2=D1*D2;
	int pt;

	for (int i=1; i<=NP; i++){
		//reshape linear array into 3x3 array: i -> (P1,P2,P3).
		P[1]=(i-1)%D1+1;
		P[2]=((i-1)%D1D2)/D1 +1;
		P[3]=(i-1)/D1D2 +1;

		pt=0;
		for (int p=1; p<=3; p++)
			for (int d=-1; d<=1; d+=2){
				pt+=1;
				P2[p]=P[p]+d;
				P2[p%3+1]=P[p%3+1];
				P2[(p+1)%3+1]=P[(p+1)%3+1];

				if (P2[p]==0 || P2[p]>Ds[p]) nn[i][pt]=0;
				else {
					nn[i][pt]=P2[1]+(P2[2]-1)*D1+(P2[3]-1)*D1D2;
			}
		}
	}
}

void interp_nn(int NP, int D1, int D2, int D3, double *values, double **allvalues, int all6, int **nn0){
	//if all6==1 gives all 6 values; if all6==0, only maximum and minimum (also including central point).

	int **nn;

	if (!nn0) {
		nn=imatrix(1,NP,1,6);
		nearest_neighbours(NP,D1,D2,D3,nn);
	}
	else nn=nn0;

	#pragma omp parallel for
	for (int i=1; i<=NP; i++){
		if (all6==0){
			allvalues[i][1]=values[i];
			allvalues[i][2]=values[i];
		}
		for (int v=1; v<=6; v++){
			if (nn[i][v]==0) continue;
			else {
				if (all6==0){
					allvalues[i][1]=fmin(allvalues[i][1],0.5*(values[nn[i][v]]+values[i]));
					allvalues[i][2]=fmax(allvalues[i][2],0.5*(values[nn[i][v]]+values[i]));
				}
				else allvalues[i][v]=0.5*(values[nn[i][v]]+values[i]);
			}
		}
	}

	if (nn0== (int **)0) free_imatrix(nn,1,NP,1,6);
}

double * duplicate_dvector(double *v, long nrl, long nrh){
	double * vnew=dvector(nrl, nrh);
	for (int r=nrl; r<=nrh; r++){
		vnew[r]=v[r];
	}
	return vnew;
}

double min_v(double *v, int N){

	double min=1/0.0;
	for (int i=0; i<N; i++) if (v[i]<min) min=v[i];
	return min;
}

double max_v(double *v, int N){
//v[0...N-1];

	double max=-1/0.0;
	for (int i=0; i<N; i++) if (v[i]>max) max=v[i];
	return max;
}

#undef NSTACK
#undef SWAP
