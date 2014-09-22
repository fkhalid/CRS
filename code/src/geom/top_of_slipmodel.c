/* Finds the top of a slip model: this valuse can later be used to impose a surface-cutting thrust (vs. an blind thrust).
 */

#include <math.h>

#include "../defines.h"


void top_of_slipmodel(struct eqkfm* eqkfm, int N){

	double z0, dz_dip;	//depth of reference point, distance along dip of closest patch;
	double len_dip;	//length of patch along dip;
	double dz;	//depth of top of shallowest patch w.r.t. reference point.
	int np;

	z0=1e30;
	for (int i=0; i<N; i++){
		np=eqkfm[i].np_di*eqkfm[i].np_st;
		len_dip=eqkfm[i].W/eqkfm[i].np_di;

		//find along dip distance of center of shallowest patch:
		dz_dip=1e30;
		for (int j=1; j<=np; j++) dz_dip=fmin(dz_dip, eqkfm[i].pos_d[j]);
		if (eqkfm[i].whichfm!=2) dz=(dz_dip-0.5*len_dip)*sin(DEG2RAD*eqkfm[i].dip1);
		else dz=(dz_dip-0.5*len_dip)*sin(DEG2RAD*eqkfm[i].dip2);

		z0=fmin(z0, eqkfm[i].depth+dz);
	}

	for (int i=0; i<N; i++) eqkfm[i].top=z0;

	return;
}
