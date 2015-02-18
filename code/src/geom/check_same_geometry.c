/*
 * check_same_geometry.c
 *
 *  Created on: Feb 18, 2015
 *      Author: camcat
 */

#include "../defines.h"

int check_same_geometry(struct eqkfm *eqkfm1, int N1, struct eqkfm *eqkfm2, int N2){
/*
 * Compares slip model structures eqkfm1, eqkfm2 to see if they have same geometry.
 *
 * Input:
 *  eqkfm1, eqkfm2: slip models
 *  N1, N2: size. eqkfm1[0...N1-1], eqkfm2[0...N2-1]
 *
 * Result:
 *  Boolean value: 1 if they have the same geometry, 0 otherwise.*
 */

	if(N1!=N2) return 0;

	for (int i=0; i<N1; i++){

		//check all things that should be equal, and if they are not return 0;
		if (eqkfm1[i].np_di!=eqkfm2[i].np_di) return 0;
		if (eqkfm1[i].np_st!=eqkfm2[i].np_st) return 0;
		if (eqkfm1[i].lat!=eqkfm2[i].lat) return 0;
		if (eqkfm1[i].lon!=eqkfm2[i].lon) return 0;
		if (eqkfm1[i].depth!=eqkfm2[i].depth) return 0;
		if (eqkfm1[i].x!=eqkfm2[i].x) return 0;
		if (eqkfm1[i].y!=eqkfm2[i].y) return 0;
		if (eqkfm1[i].L!=eqkfm2[i].L) return 0;
		if (eqkfm1[i].W!=eqkfm2[i].W) return 0;

		for (int j=1; j<=eqkfm1[i].np_di*eqkfm1[i].np_st; j++){
			if (eqkfm1[i].pos_d[j]!=eqkfm2[i].pos_d[j]) return 0;
			if (eqkfm1[i].pos_s[j]!=eqkfm2[i].pos_s[j]) return 0;
		}

	}

	//everything is equal, hence return 1;
	return 1;
}
