# include "ran1.h"

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
