/*
 * GR.c
 *
 *  Created on: Sep 20, 2013
 *      Author: camcat
 */

#include "GR.h"

double *assign_GRnorm(double *mags, int N, double b, int Minf){
/* mags are vector containing the center of sorted, equally spaced magnitude bins. Indices: [1...N].
 * GRpar={a,b} of GR distrib: N(M>=m)=10^(a-b*m).
 * Minf is flag indicating if last bin is open *
 */

	double dm;
	double norm;
	double *weights;
	double one;
	double mmax, mmin;
	double a=1;	//doesn't matter since weights are normalized.

	if (N==1) {
		weights=dvector(1,1);
		weights[1]=1;
	}

	else {
		weights=dvector(1,N);
		dm=mags[2]-mags[1];
		mmin=mags[1]-0.5*dm;
		mmax=mags[N]+0.5*dm;
		norm=Minf? pow(10,a-b*mmin) : pow(10,a-b*mmin)-pow(10,a-b*mmax);
		for (int i=1; i<=N; i++){
			mmin=mags[i]-0.5*dm;
			mmax=mags[i]+0.5*dm;
			weights[i]=(pow(10,a-b*mmin)-pow(10,a-b*mmax))/norm;
		}
		if (Minf) weights[N]=(pow(10,a-b*mmin))/norm;
	}

	return weights;

}

double Mc_maxcurv(double *mags, int N){
//find magnitude of completeness using maximum curvature method; events are
//binned so that each bin has same no of counts, and normalized by bin size.
//vector mags has indices [0...N-1].

	double *bin_c, *bin_h;
	double f=0.0;
	double *mags2;
	int Nbin=floor(sqrt(N)), ind;

	mags2=dvector(0,N-1);
	for (int i=0; i<N; i++) mags2[i]=mags[i];	//since next function sorts array.

	bin_equnumber(mags2,N, Nbin, &bin_c, &bin_h);

	for (int i=0; i<Nbin; i++){
		if (f<=bin_h[i]){
			f=bin_h[i];
			ind=i;
		}
	}

	free_dvector(mags2, 0, N);
	return bin_c[ind];
}

double calculatebvalue(double *mags, int N, double Mc){
  double  *mags2;
  double dM, dum, mean, b, db;
  int Z;

  mags2=dvector(0,N-1);

  Z=0;
  mean=0.0;
  //pick magnitudes above completeness, calculate mean magnitude and fill in new vector.
  for (int i=0; i<N; i++){
	  if (mags[i]>=Mc){
		  mags2[Z]=mags[i];
		  mean+=mags[i];
		  Z+=1;
	  }
  }
  mean*=(1.0/Z);

  qsort (mags2, Z, sizeof(double), &compare);

  dM=1000.0;
  for(int i=1;i<Z;i++) {
	  dum=(mags2[i]-mags2[i-1]);
	  if(dum!=0.0 && dum<dM) dM=dum;
  }

// Maximum Likelihood fit result:
  b  = log10(exp(1.0))/(mean-(Mc-0.5*dM));
  db = 1.96*b/sqrt(1.0*Z);

  return b;
}

int bin_equnumber(double *v, int N, int Nbin, double **bin_c, double **norm_count){
//norm_count is the no of events per bin divided by bin size.

	double bot, top;
	int Nel, Nellast;

	qsort (v, N, sizeof(double), &compare);
	Nel=floor(N/Nbin);

	*bin_c=dvector(0,Nbin-1);
	*norm_count=dvector(0,Nbin-1);


	for (int i=0; i<Nbin-1; i++){

		bot= (i==0) ? v[0] : 0.5*(v[Nel*i]+v[Nel*i-1]);
		top= 0.5*(v[Nel*(i+1)-1]+v[Nel*(i+1)]);

		(*bin_c)[i]=0.5*(top+bot);
		(*norm_count)[i]=Nel/(top-bot);
	}

	bot=0.5*(v[Nel*(Nbin-1)]+v[Nel*(Nbin-1)-1]);
	top=v[N-1];
	Nellast=N-Nel*(Nbin-1);
	(*bin_c)[Nbin-1]=0.5*(top+bot);
	(*norm_count)[Nbin-1]=Nellast/(top-bot);

	return (0);

}

int compare (const void * a, const void * b)
{
	if (*(double*)a==*(double*)b) return 0;
	else if (*(double*)a>*(double*)b) return 1;
	else return -1;
}
