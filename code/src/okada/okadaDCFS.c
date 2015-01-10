/*
 * okadaDCFS.c
 *
 *  Created on: Jan 4, 2013
 *      Author: camcat
 */

#include "okadaDCFS.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

//todo make all these functions return error code.
//todo check if there are memory leaks.

//---------------------------------------------------------------------//
//-----					Top level functions						  -----//
//---------------------------------------------------------------------//

int resolve_DCFS(struct pscmp DCFS, struct crust crst, double *strikeRs, double *dipRs, double *rake, int optrake){
	//To be called after filling in DCFS.S. (with functions below).
	//todo make vectors static to avoid repeated memory allocations?

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double *sigma0s;
	double **stress0s;
	double **n, **s;
	double MaxDCFS=DCFS_cap;
	int Nsel=DCFS.nsel;
	int fm;
	double strikeR, dipR;
	int no_fm_zones=crst.nofmzones;

	n=malloc((no_fm_zones)*sizeof(double *));
	s=malloc((no_fm_zones)*sizeof(double *));
	stress0s=malloc((no_fm_zones)*sizeof(double *));
	sigma0s=dvector(0,no_fm_zones-1);

	for (int i=0; i<no_fm_zones; i++) {
		strikeR=strikeRs[i];
		dipR=dipRs[i];
		n[i]=normal_vector(strikeR, dipR);

		if (!optrake){
			if (!rake) {
				print_screen("** Warning: optrake=0, but rake is NULL: will use optimal rake (resolve_DCFS).**\n");
				print_logfile("** Warning: optrake=0, but rake is NULL: will use optimal rake (resolve_DCFS).**\n");
				optrake=1;
			}
			s[i]=slip_vector(strikeR, dipR, *rake);
		}
		else s[i]=NULL;

		//calculated background stress and normal stress:
		stress0s[i]=mtimesv(crst.S,n[i],NULL,3,3);
		sigma0s[i]=-1.0*vdotv(stress0s[i],n[i],3);		// compression is positive (rock mechanics convention).
	}

	#pragma omp parallel for private(fm)
	for (int i=1; i<=Nsel; i++){
		fm= (no_fm_zones==1) ? 0 : crst.fmzone[i];
		DCFS.cmb[i]=resolve_n(DCFS.S[i], n[fm], NULL, crst.fric, stress0s[fm], sigma0s[fm], s[fm]);
		if (DCFS.cmb[i]>MaxDCFS) DCFS.cmb[i]=MaxDCFS;
		if (DCFS.cmb[i]<-MaxDCFS) DCFS.cmb[i]=-MaxDCFS;
	}

	for (int i=0; i<no_fm_zones; i++){
		if (s[i]) free_dvector(s[i],1,3);
		free_dvector(n[i],1,3);
		free_dvector(stress0s[i],1,3);
	}
	free(s);
	free(n);
	free(stress0s);
	free_dvector(sigma0s,0,1);

	return(0);
}

// todo [coverage] this block is never tested
int okadaDCFS(struct pscmp DCFS, struct eqkfm *eqkfm1, int NF, struct crust crst, double *strikeR, double *dipR, int full_tensor){

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double north, east, eqnorth, eqeast;
	double len, width, depth; //for individual patches.
	double strike, dip, rake;
	double alpha;
	double Sxx, Syy, Szz, Sxy, Syz, Sxz;
	int i2, Nsel, err=0;

	alpha = (crst.lambda + crst.mu)/(crst.lambda + 2*crst.mu);

	if (DCFS.nsel!=(*eqkfm1).nsel){
		print_screen("**Warning: DCFS.nsel!=eqkfm.nsel in okadaDCFSc.**\n");
		print_logfile("**Error: DCFS.nsel!=eqkfm.nsel in okadaDCFSc. Will choose the one with more points selected.**\n");
		DCFS.which_pts=(*eqkfm1).selpoints= (DCFS.nsel>(*eqkfm1).nsel) ? DCFS.which_pts : (*eqkfm1).selpoints;
		DCFS.nsel=(*eqkfm1).nsel=fmax(DCFS.nsel,(*eqkfm1).nsel);
	}

	Nsel=DCFS.nsel;

	//---------initialize DCFS----------//

	for (int i=1; i<=DCFS.nsel; i++){
		DCFS.cmb[i]=0;
		for (int a=1; a<=3; a++){
			for (int b=1; b<=3; b++) DCFS.S[i][a][b]=0;
		}
	}

	//-----------------------------------------------------------------------------------------//
	//-----------Calculate Coulomb stress vector for each patch and add them up.---------------//
	//-----------------------------------------------------------------------------------------//

	for (int j=0; j<NF; j++){

		//patch size:
		len=eqkfm1[j].L*(1.0/eqkfm1[j].np_st);
		width=eqkfm1[j].W*(1.0/eqkfm1[j].np_di);

		if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
			print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);

			return(1);
		}

		for (int p=1; p<=eqkfm1[j].np_di*eqkfm1[j].np_st; p++) {
			// calculate patch position:
			patch_pos(eqkfm1[j], p, &eqeast, &eqnorth, &depth);

			// calculate DCFS from patches:
			//#pragma omp parallel for private(Sxx, Syy, Szz, Sxy, Syz, Sxz, north, east, i2)
			for (int i=1; i<=Nsel; i++){
				i2=DCFS.which_pts[i];
				north=crst.y[i2];
				east=crst.x[i2];
				pscokada(eqnorth, eqeast, depth, strike, dip, len, width, eqkfm1[j].slip_str[p],
						 -1.0*eqkfm1[j].slip_dip[p], north, east, crst.depth[i2], &Sxx, &Syy,
						 &Szz, &Sxy, &Syz, &Sxz, alpha, crst.lambda, crst.mu, crst.fric);
				DCFS.S[i][1][1]+=1e6*Sxx;
				DCFS.S[i][2][2]+=1e6*Syy;
				DCFS.S[i][3][3]+=1e6*Szz;
				DCFS.S[i][1][2]+=1e6*Sxy;
				DCFS.S[i][2][3]+=1e6*Syz;
				DCFS.S[i][1][3]+=1e6*Sxz;
			}
		}
	}

	for (int i=1; i<=Nsel; i++) {
		DCFS.S[i][2][1]=DCFS.S[i][1][2];
		DCFS.S[i][3][2]=DCFS.S[i][2][3];
		DCFS.S[i][3][1]=DCFS.S[i][1][3];
	}

	//-----------------------------------------------------------------------------------------//
	//-----------------Fill in Syx, Szy, Szx and resolve stress on required plane.-------------//
	//-----------------------------------------------------------------------------------------//

	if (!full_tensor) resolve_DCFS(DCFS, crst, strikeR, dipR, NULL, 1);	//todo not hardwire optrake?

	return(0);
}

//int okadaCoeff_mpi(float ****Coeffs_st, float ****Coeffs_dip, struct eqkfm *eqkfm1,
//			   int NF, struct crust crst, double *lats, double *lons, double *depths) {
//	//lats, lons, depths contain complete list of grid points.
//	// Only the ones with indices eqkfm1.selpoints will be used.
//
//	// [Fahad] Variables used for MPI.
//	int procId = 0, numProcs = 1;
//	int start, end, partitionSize;
//
//	#ifdef _CRS_MPI
//		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
//		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
//	#endif
//
//	double north, east, eqnorth, eqeast;
//	double len, width, depth; //for individual patches.
//	double depth0; //to differentiate between blind fault (depth0=0) or fault cutting through surface.
//	double strike, dip, rake;
//	double alpha;
//	double Sxx, Syy, Szz, Sxy, Syz, Sxz;
//	int Nsel = eqkfm1[0].nsel;
//	int NP_tot = 0, p1, i;
//	int pure_thrustnorm, pure_strslip;
//	int err=0;
//
//	for(int j=0; j<NF; j++) {
//		NP_tot+=eqkfm1[j].np_di*eqkfm1[j].np_st;
//	}
//
//	alpha = (crst.lambda + crst.mu)/(crst.lambda + 2*crst.mu);
//	depth0=eqkfm1[0].cuts_surf ? eqkfm1[0].top : 0.0;
//
//	print_logfile("Depth of surface: %.3lf km.\n", depth0);
//
//	//---------initialize DCFS----------//
//	*Coeffs_st  = f3tensor(1, NP_tot, 1, Nsel, 1, 6);	//TODO should deallocate at the end (in main.c).
//	*Coeffs_dip = f3tensor(1, NP_tot, 1, Nsel, 1, 6);
////
////	for (int p1=1; p1<=NP_tot-1; p1++){
////		for (int i=1; i<=Nsel; i++){
////			for (int j=1; j<=6; j++){
////				(*Coeffs_st)[p1][i][j]=0.0;
////				(*Coeffs_dip)[p1][i][j]=0.0;
////			}
////		}
////	}
////	float ***Coeffs_st1, ***Coeffs_dip1;
////	Coeffs_st1  = f3tensor(1,NP_tot,1,Nsel,1,6);
////	Coeffs_dip1 = f3tensor(1,NP_tot,1,Nsel,1,6);
//
////	for (int p1=1; p1<=NP_tot; p1++){
////		for (int i=1; i<=Nsel; i++){
////			for (int j=1; j<=6; j++){
////				(*Coeffs_st)[p1][i][j]=0.0;
////				(*Coeffs_dip)[p1][i][j]=0.0;
////			}
////		}
////	}
//
//
//	// [Fahad] - Create a linearized full tensor.
////	size_t fullTensorSize = ((NP_tot-1) * Nsel * 6) - 1;
//	size_t fullTensorSize = ((NP_tot) * Nsel * 6);
//	float *coeffs_st  = (float*) malloc((size_t)(fullTensorSize * sizeof(float)));
//	float *coeffs_dip = (float*) malloc((size_t)(fullTensorSize * sizeof(float)));
//
//	//-----------------------------------------------------------------------------------------//
//	//-----------Calculate Coulomb stress vector for each patch assuming slip=1.---------------//
//	//-----------------------------------------------------------------------------------------//
//
//	print_screen("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);
//	print_logfile("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);
//
//	p1=0;	//count total number of patches (in all faults);
//	for (int j=0; j<NF; j++) {
//		pure_thrustnorm=pure_strslip=0;
//
//		if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
//			print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
//			print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
//			return(1);
//		}
//
//		len=eqkfm1[j].L*(1.0/eqkfm1[j].np_st);
//		width=eqkfm1[j].W*(1.0/eqkfm1[j].np_di);
//
////		int numPatches = eqkfm1[j].np_di*eqkfm1[j].np_st - 1;
//		int numPatches = eqkfm1[j].np_di*eqkfm1[j].np_st;
//
//		#ifdef _CRS_MPI
//			// FIXME: Write a simple algorithm to fit lower Nsur values to numProcs ...
//			if(numProcs > numPatches) {
//				if(procId == 0) {
//					printf("\n Number of processes: %d", numProcs);
//					printf("\n Number of iterations: %d", numPatches);
//				}
//				error_quit("\n **numPatches must be greater than or equal to the number of processes ** \n");
//			}
//
//			partitionSize = numPatches / numProcs;
//
//			// If partionSize leaves one or more elements at the end
//			if(partitionSize * numProcs != numPatches) {
//				partitionSize += 1;
//
//				coeffs_st  = (float*) realloc(coeffs_st,  (size_t)((partitionSize*Nsel*6*numProcs) * sizeof(float)));
//				coeffs_dip = (float*) realloc(coeffs_dip, (size_t)((partitionSize*Nsel*6*numProcs) * sizeof(float)));
//			}
//
//			start = (procId * partitionSize);
//			end = start + partitionSize;
//
////			// The very last process
////			if(procId / (numProcs-1) == 1) {
////				// If partionSize leaves one or more elements at the end
////				if(partitionSize * numProcs != numPatches) {
////					// Include the remaining elements
////					partitionSize += numPatches - (partitionSize * numProcs);
////				}
////			}
//
////			if(procId == 0) {
////				start = 1;
////			}
//		#else
//			start = 1;
//			end = numPatches + 1;
//		#endif
//
//		// [Fahad] - Create a linearized partitioned tensor.
//		size_t partitionedTensorSize = partitionSize * Nsel * 6;
//		float *coeffs_st_partitioned  = (float*) malloc((size_t)(partitionedTensorSize * sizeof(float)));
//		float *coeffs_dip_partitioned = (float*) malloc((size_t)(partitionedTensorSize * sizeof(float)));
//
//		for(int i = 0; i < partitionedTensorSize; ++i) {
//			coeffs_st_partitioned [i] = 0.0;
//			coeffs_dip_partitioned[i] = 0.0;
//		}
//
////		if(procId != 0) {
////			printf("\n numPatches: %d", numPatches);
////			printf("\n partitionSize: %d", partitionSize);
////			printf("\n partitionedTensorSize: %d", partitionedTensorSize);
////		}
//
//		int p2 = start, index=0;
//		for(int p = 0; p < partitionSize; p++) {
//			patch_pos(eqkfm1[j], p2+1, &eqeast, &eqnorth, &depth);
//			++p2;
//
//			#pragma omp parallel for private(Sxx, Syy, Szz, Sxy, Syz, Sxz, north, east, i)
//			for(int i0=0; i0<Nsel; i0++) {
//				i=eqkfm1[0].selpoints[i0+1];	// [Fahad] Added '1' to the index
//				north=crst.y[i];
//				east=crst.x[i];
//
//				if(pure_thrustnorm!=1) {
//					pscokada(eqnorth, eqeast, depth-depth0,  strike,  dip, len, width, 1, 0,
//							north, east, depths[i]-depth0, &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
//							alpha, crst.lambda, crst.mu, crst.fric);
//
//					index = (p * Nsel * 6) + (i0 * 6);
//
//					coeffs_st_partitioned[index + 0] += 1e6*Sxx;
//					coeffs_st_partitioned[index + 1] += 1e6*Syy;
//					coeffs_st_partitioned[index + 2] += 1e6*Szz;
//					coeffs_st_partitioned[index + 3] += 1e6*Sxy;
//					coeffs_st_partitioned[index + 4] += 1e6*Syz;
//					coeffs_st_partitioned[index + 5] += 1e6*Sxz;
//				}
//
//				if(pure_strslip!=1) {
//					pscokada(eqnorth, eqeast, depth-depth0,  strike, dip, len, width, 0, -1,
//							 north, east, depths[i]-depth0, &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
//							 alpha, crst.lambda, crst.mu, crst.fric);
//
//					index = (p * Nsel * 6) + (i0 * 6);
//
//					coeffs_dip_partitioned[index + 0] += 1e6*Sxx;
//					coeffs_dip_partitioned[index + 1] += 1e6*Syy;
//					coeffs_dip_partitioned[index + 2] += 1e6*Szz;
//					coeffs_dip_partitioned[index + 3] += 1e6*Sxy;
//					coeffs_dip_partitioned[index + 4] += 1e6*Syz;
//					coeffs_dip_partitioned[index + 5] += 1e6*Sxz;
//				}
//			}
//		}
//
////		if(procId != 0) {
////			printf("\n ProcId %d -- coeffs_st_partitioned[0]:  %f \n", procId, coeffs_st_partitioned[0]);
////			printf("\n ProcId %d -- coeffs_st_partitioned[(partitionSize-1)*Nsel*6]:  %f \n", procId, coeffs_st_partitioned[(partitionSize-1)*Nsel*6]);
////		}
//
//		#ifdef _CRS_MPI
//			MPI_Allgather(coeffs_st_partitioned, partitionedTensorSize,
//						  MPI_FLOAT, coeffs_st, partitionedTensorSize,
//						  MPI_FLOAT, MPI_COMM_WORLD);
//
//			MPI_Allgather(coeffs_dip_partitioned, partitionedTensorSize,
//						  MPI_FLOAT, coeffs_dip, partitionedTensorSize,
//						  MPI_FLOAT, MPI_COMM_WORLD);
//		#endif
//
//		free(coeffs_st_partitioned);
//		free(coeffs_dip_partitioned);
//	}
//
////	if(procId != 0) {
////		printf("\n ProcId %d -- coeffs_st[0]:  %f \n", procId, coeffs_st[0]);
////		printf("\n ProcId %d -- coeffs_st[80*Nsel*6]:  %f \n", procId, coeffs_st[80*Nsel*6]);
////		printf("\n ProcId %d -- coeffs_st[81*Nsel*6]:  %f \n", procId, coeffs_st[81*Nsel*6]);
////		printf("\n ProcId %d -- coeffs_st[160*Nsel*6]:  %f \n", procId, coeffs_st[160*Nsel*6]);
////	}
//
//	// [Fahad] - Copy data from the linearized tensor to the f3tensor.
//	int index = 0;
//	for(int i = 0; i < NP_tot; ++i) {
//		for(int j = 0; j < Nsel; ++j) {
//			for(int k = 0; k < 6; ++k) {
//				index = (i * Nsel * 6) + (j * 6) + k;
//
//				(*Coeffs_st) [i+1][j+1][k+1] = coeffs_st [index];
//				(*Coeffs_dip)[i+1][j+1][k+1] = coeffs_dip[index];
//			}
//		}
//	}
//
//	if(procId != 0) {
////		printf("\n\n ProcId %d -- (*Coeffs_st) [1][1][1]:  %f \n", procId, (*Coeffs_st1)[1][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [81][1][1]: %f \n", procId, (*Coeffs_st1)[81][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [82][1][1]: %f \n", procId, (*Coeffs_st1)[82][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [83][1][1]: %f \n", procId, (*Coeffs_st1)[83][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [161][1][1]: %f \n", procId, (*Coeffs_st1)[161][1][1]);
////
////		printf("\n\n ProcId %d -- (*Coeffs_dip) [1][1][1]:  %f \n", procId, (*Coeffs_dip)[1][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip) [81][1][1]: %f \n", procId, (*Coeffs_dip)[81][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip) [82][1][1]: %f \n", procId, (*Coeffs_dip)[82][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip) [83][1][1]: %f \n", procId, (*Coeffs_dip)[83][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip) [161][1][1]: %f \n", procId, (*Coeffs_dip)[161][1][6]);
//
////		printf("\n ProcId %d -- (*Coeffs_st) [1][1][1]:  %f \n", procId,  (*Coeffs_st)[1][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [94][1][1]: %f \n", procId,  (*Coeffs_st)[94][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [95][1][1]: %f \n", procId,  (*Coeffs_st)[95][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [96][1][1]: %f \n", procId,  (*Coeffs_st)[96][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_st) [188][1][1]: %f \n", procId, (*Coeffs_st)[188][1][1]);
////
////		printf("\n\n ProcId %d -- (*Coeffs_dip)[1][1][1]:  %f \n", procId, (*Coeffs_dip)[1][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip)[94][1][1]: %f \n", procId,   (*Coeffs_dip)[94][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip)[95][1][1]: %f \n", procId,   (*Coeffs_dip)[95][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip)[96][1][1]: %f \n", procId,   (*Coeffs_dip)[96][1][1]);
////		printf("\n ProcId %d -- (*Coeffs_dip)[188][1][1]: %f \n", procId,  (*Coeffs_dip)[188][1][1]);
//	}
//
//	free(coeffs_st);
//	free(coeffs_dip);
//
////	for (int p1=1; p1<=NP_tot; p1++){
////		for (int i=1; i<=Nsel; i++){
////			for (int j=1; j<=6; j++){
////				if((*Coeffs_st)[p1][i][j] != Coeffs_st1[p1][i][j]) {
////					if(procId != 0) {
////						printf("\n Mismatch detected: %f --- %f", (*Coeffs_st)[p1][i][j], Coeffs_st1[p1][i][j]);
////					}
////				}
////			}
////		}
////	}
////
////	(*Coeffs_st) = Coeffs_st1;
//
//	return(0);
//}

// Same function as the one commented out above, just cleaned up a bit.
int okadaCoeff_mpi(float ****Coeffs_st, float ****Coeffs_dip, struct eqkfm *eqkfm1,
			   int NF, struct crust crst, double *lats, double *lons, double *depths) {
	//lats, lons, depths contain complete list of grid points.
	// Only the ones with indices eqkfm1.selpoints will be used.

	// [Fahad] Variables used for MPI.
	int procId = 0, numProcs = 1;
	int start, end, partitionSize;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	#endif

	double north, east, eqnorth, eqeast;
	double len, width, depth; //for individual patches.
	double depth0; //to differentiate between blind fault (depth0=0) or fault cutting through surface.
	double strike, dip, rake;
	double alpha;
	double Sxx, Syy, Szz, Sxy, Syz, Sxz;
	int Nsel = eqkfm1[0].nsel;
	int NP_tot = 0, p1, i;
	int pure_thrustnorm, pure_strslip;
	int err=0;

	for(int j=0; j<NF; j++) {
		NP_tot+=eqkfm1[j].np_di*eqkfm1[j].np_st;
	}

	alpha = (crst.lambda + crst.mu)/(crst.lambda + 2*crst.mu);
	depth0=eqkfm1[0].cuts_surf ? eqkfm1[0].top : 0.0;

	print_logfile("Depth of surface: %.3lf km.\n", depth0);

	//---------initialize DCFS----------//
	*Coeffs_st  = f3tensor(1, NP_tot, 1, Nsel, 1, 6);	//TODO should deallocate at the end (in main.c).
	*Coeffs_dip = f3tensor(1, NP_tot, 1, Nsel, 1, 6);

	// [Fahad] - Create a linearized full tensor.
	size_t fullTensorSize = ((NP_tot) * Nsel * 6);
	float *coeffs_st  = (float*) malloc((size_t)(fullTensorSize * sizeof(float)));
	float *coeffs_dip = (float*) malloc((size_t)(fullTensorSize * sizeof(float)));

	//-----------------------------------------------------------------------------------------//
	//-----------Calculate Coulomb stress vector for each patch assuming slip=1.---------------//
	//-----------------------------------------------------------------------------------------//

	print_screen("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);
	print_logfile("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);

	p1=0;	//count total number of patches (in all faults);
	for (int j=0; j<NF; j++) {
		pure_thrustnorm=pure_strslip=0;

		if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
			print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			return(1);
		}

		len=eqkfm1[j].L*(1.0/eqkfm1[j].np_st);
		width=eqkfm1[j].W*(1.0/eqkfm1[j].np_di);

		int numPatches = eqkfm1[j].np_di*eqkfm1[j].np_st;

		#ifdef _CRS_MPI
			if(numProcs > numPatches) {
				if(procId == 0) {
					printf("\n Number of processes: %d", numProcs);
					printf("\n Number of iterations: %d", numPatches);
				}
				error_quit("\n **numPatches must be greater than or equal to the number of processes ** \n");
			}

			partitionSize = numPatches / numProcs;

			// If partionSize leaves one or more elements at the end
			if(partitionSize * numProcs != numPatches) {
				partitionSize += 1;

				coeffs_st  = (float*) realloc(coeffs_st,  (size_t)((partitionSize*Nsel*6*numProcs) * sizeof(float)));
				coeffs_dip = (float*) realloc(coeffs_dip, (size_t)((partitionSize*Nsel*6*numProcs) * sizeof(float)));
			}

			start = (procId * partitionSize);
			end = start + partitionSize;
		#else
			start = 1;
			end = numPatches + 1;
		#endif

		// [Fahad] - Create a linearized partitioned tensor.
		size_t partitionedTensorSize = partitionSize * Nsel * 6;
		float *coeffs_st_partitioned  = (float*) malloc((size_t)(partitionedTensorSize * sizeof(float)));
		float *coeffs_dip_partitioned = (float*) malloc((size_t)(partitionedTensorSize * sizeof(float)));

		for(int i = 0; i < partitionedTensorSize; ++i) {
			coeffs_st_partitioned [i] = 0.0;
			coeffs_dip_partitioned[i] = 0.0;
		}

		int p2 = start, index=0;
		for(int p = 0; p < partitionSize; p++) {
			patch_pos(eqkfm1[j], p2+1, &eqeast, &eqnorth, &depth);
			++p2;

			#pragma omp parallel for private(Sxx, Syy, Szz, Sxy, Syz, Sxz, north, east, i)
			for(int i0=0; i0<Nsel; i0++) {
				i=eqkfm1[0].selpoints[i0+1];	// [Fahad] Added '1' to the index
				north=crst.y[i];
				east=crst.x[i];

				if(pure_thrustnorm!=1) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike,  dip, len, width, 1, 0,
							north, east, depths[i]-depth0, &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
							alpha, crst.lambda, crst.mu, crst.fric);

					index = (p * Nsel * 6) + (i0 * 6);

					coeffs_st_partitioned[index + 0] += 1e6*Sxx;
					coeffs_st_partitioned[index + 1] += 1e6*Syy;
					coeffs_st_partitioned[index + 2] += 1e6*Szz;
					coeffs_st_partitioned[index + 3] += 1e6*Sxy;
					coeffs_st_partitioned[index + 4] += 1e6*Syz;
					coeffs_st_partitioned[index + 5] += 1e6*Sxz;
				}

				if(pure_strslip!=1) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike, dip, len, width, 0, -1,
							 north, east, depths[i]-depth0, &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
							 alpha, crst.lambda, crst.mu, crst.fric);

					index = (p * Nsel * 6) + (i0 * 6);

					coeffs_dip_partitioned[index + 0] += 1e6*Sxx;
					coeffs_dip_partitioned[index + 1] += 1e6*Syy;
					coeffs_dip_partitioned[index + 2] += 1e6*Szz;
					coeffs_dip_partitioned[index + 3] += 1e6*Sxy;
					coeffs_dip_partitioned[index + 4] += 1e6*Syz;
					coeffs_dip_partitioned[index + 5] += 1e6*Sxz;
				}
			}
		}

		#ifdef _CRS_MPI
			MPI_Allgather(coeffs_st_partitioned, partitionedTensorSize,
						  MPI_FLOAT, coeffs_st, partitionedTensorSize,
						  MPI_FLOAT, MPI_COMM_WORLD);

			MPI_Allgather(coeffs_dip_partitioned, partitionedTensorSize,
						  MPI_FLOAT, coeffs_dip, partitionedTensorSize,
						  MPI_FLOAT, MPI_COMM_WORLD);
		#endif

		free(coeffs_st_partitioned);
		free(coeffs_dip_partitioned);
	}

	// [Fahad] - Copy data from the linearized tensor to the f3tensor.
	int index = 0;
	for(int i = 0; i < NP_tot; ++i) {
		for(int j = 0; j < Nsel; ++j) {
			for(int k = 0; k < 6; ++k) {
				index = (i * Nsel * 6) + (j * 6) + k;

				(*Coeffs_st) [i+1][j+1][k+1] = coeffs_st [index];
				(*Coeffs_dip)[i+1][j+1][k+1] = coeffs_dip[index];
			}
		}
	}

	free(coeffs_st);
	free(coeffs_dip);

	return(0);
}

// todo [coverage] this block is never tested
int okadaCoeff(float ****Coeffs_st, float ****Coeffs_dip, struct eqkfm *eqkfm1, int NF,
			   struct crust crst, double *lats, double *lons, double *depths) {
	//lats, lons, depths contain complete list of grid points. Only the ones with indices eqkfm1.selpoints will be used.

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double north, east, eqnorth, eqeast;
	double len, width, depth; //for individual patches.
	double strike, dip, rake;
	double alpha;
	double Sxx, Syy, Szz, Sxy, Syz, Sxz;
	int Nsel=eqkfm1[0].nsel;
	int NP_tot=0, p1, i;
	int pure_thrustnorm, pure_strslip;
	int err=0;

	for (int j=0; j<NF; j++) NP_tot+=eqkfm1[j].np_di*eqkfm1[j].np_st;

	alpha = (crst.lambda + crst.mu)/(crst.lambda + 2*crst.mu);

	//---------initialize DCFS----------//

	*Coeffs_st=f3tensor(1,NP_tot,1,Nsel,1,6);	//TODO should deallocate at the end (in main.c).
	*Coeffs_dip=f3tensor(1,NP_tot,1,Nsel,1,6);

	for (int p1=1; p1<=NP_tot; p1++){
		for (int i=1; i<=Nsel; i++){
			for (int j=1; j<=6; j++){
				(*Coeffs_st)[p1][i][j]=0.0;
				(*Coeffs_dip)[p1][i][j]=0.0;
			}
		}
	}

	//-----------------------------------------------------------------------------------------//
	//-----------Calculate Coulomb stress vector for each patch assuming slip=1.---------------//
	//-----------------------------------------------------------------------------------------//

	print_screen("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);
	print_logfile("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);

	p1=0;	//count total number of patches (in all faults);
	for (int j=0; j<NF; j++){
		pure_thrustnorm=pure_strslip=0;
//		todo: check rake for all patches:
//		if (fmod(eqkfm1[j].rake1+90,180.0)!=0) pure_thrustnorm=0;
//		if (fmod(eqkfm1[j].rake1,180.0)!=0) pure_strslip=0;

		if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
			print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			return(1);
		}

		len=eqkfm1[j].L*(1.0/eqkfm1[j].np_st);
		width=eqkfm1[j].W*(1.0/eqkfm1[j].np_di);

		for (int p=1; p<=eqkfm1[j].np_di*eqkfm1[j].np_st; p++){
			p1+=1;
			patch_pos(eqkfm1[j], p, &eqeast, &eqnorth, &depth);

			#pragma omp parallel for private(Sxx, Syy, Szz, Sxy, Syz, Sxz, north, east, i)
			for (int i0=1; i0<=Nsel; i0++){
				i=eqkfm1[0].selpoints[i0];
				north=crst.y[i];
				east=crst.x[i];
				if (pure_thrustnorm!=1) {
					pscokada(eqnorth, eqeast, depth,  strike,  dip, len, width, 1, 0,
							 north, east, depths[i], &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
							 alpha, crst.lambda, crst.mu, crst.fric);

					(*Coeffs_st)[p1][i0][1]+=1e6*Sxx;
					(*Coeffs_st)[p1][i0][2]+=1e6*Syy;
					(*Coeffs_st)[p1][i0][3]+=1e6*Szz;
					(*Coeffs_st)[p1][i0][4]+=1e6*Sxy;
					(*Coeffs_st)[p1][i0][5]+=1e6*Syz;
					(*Coeffs_st)[p1][i0][6]+=1e6*Sxz;
				}

				if (pure_strslip!=1){
					pscokada(eqnorth, eqeast, depth,  strike, dip, len, width, 0, -1,
							 north, east, depths[i], &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
							 alpha, crst.lambda, crst.mu, crst.fric);

					(*Coeffs_dip)[p1][i0][1]+=1e6*Sxx;
					(*Coeffs_dip)[p1][i0][2]+=1e6*Syy;
					(*Coeffs_dip)[p1][i0][3]+=1e6*Szz;
					(*Coeffs_dip)[p1][i0][4]+=1e6*Sxy;
					(*Coeffs_dip)[p1][i0][5]+=1e6*Syz;
					(*Coeffs_dip)[p1][i0][6]+=1e6*Sxz;
				}
			}
		}

		if(procId == 0) {
//			printf("\n\n ProcId %d -- (*Coeffs_st) [1][1][1]:  %f \n", procId, (*Coeffs_st)[1][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [81][1][1]: %f \n", procId, (*Coeffs_st)[81][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [82][1][1]: %f \n", procId, (*Coeffs_st)[82][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [83][1][1]: %f \n", procId, (*Coeffs_st)[83][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [161][1][1]: %f \n", procId, (*Coeffs_st)[161][1][1]);
//
//			printf("\n\n ProcId %d -- (*Coeffs_st) [1][1][1]:  %f \n", procId, (*Coeffs_dip)[1][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [81][1][1]: %f \n", procId, (*Coeffs_dip)[81][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [82][1][1]: %f \n", procId, (*Coeffs_dip)[82][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [83][1][1]: %f \n", procId, (*Coeffs_dip)[83][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [161][1][1]: %f \n", procId, (*Coeffs_dip)[161][1][6]);

//			printf("\n ProcId %d -- (*Coeffs_st) [1][1][1]:  %f \n", procId, (*Coeffs_st)[1][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [94][1][1]: %f \n", procId, (*Coeffs_st)[94][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [95][1][1]: %f \n", procId, (*Coeffs_st)[95][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [96][1][1]: %f \n", procId, (*Coeffs_st)[96][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [188][1][1]: %f \n", procId, (*Coeffs_st)[188][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[1][1][1]:  %f \n", procId, (*Coeffs_dip)[1][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[94][1][1]: %f \n", procId, (*Coeffs_dip)[94][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[95][1][1]: %f \n", procId, (*Coeffs_dip)[95][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[96][1][1]: %f \n", procId, (*Coeffs_dip)[96][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [187][1][1]: %f \n", procId, (*Coeffs_st)[187][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [188][1][1]: %f \n", procId, (*Coeffs_st)[188][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [189][1][1]: %f \n", procId, (*Coeffs_st)[189][1][1]);


//			printf("\n ProcId %d -- (*Coeffs_st) [1][1][1]:  %f \n", procId, (*Coeffs_st)[1][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [94][1][1]: %f \n", procId, (*Coeffs_st)[94][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [95][1][1]: %f \n", procId, (*Coeffs_st)[95][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_st) [189][1][1]: %f \n", procId, (*Coeffs_st)[188][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[1][1][1]:  %f \n", procId, (*Coeffs_dip)[1][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[94][1][1]: %f \n", procId, (*Coeffs_dip)[94][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[95][1][1]: %f \n", procId, (*Coeffs_dip)[95][1][1]);
//			printf("\n ProcId %d -- (*Coeffs_dip)[189][1][1]: %f \n", procId, (*Coeffs_dip)[188][1][1]);
		}

//		printf("\n ProcId %d -- (*Coeffs_st)[1][1][1]: %f \n", procId, (*Coeffs_st)[1][1][1]);
//		printf("\n ProcId %d -- (*Coeffs_st)[280][1][1]: %f \n", procId, (*Coeffs_st)[280][1][1]);
//		printf("\n ProcId %d -- (*Coeffs_st)[281][1][1]: %f \n", procId, (*Coeffs_st)[281][1][1]);
//		printf("\n ProcId %d -- (*Coeffs_st)[560][1][1]: %f \n", procId, (*Coeffs_st)[560][1][1]);
//		printf("\n ProcId %d -- (*Coeffs_dip)[1][1][1]: %f \n", procId, (*Coeffs_dip)[1][1][1]);
//		printf("\n ProcId %d -- (*Coeffs_dip)[280][1][1]: %f \n", procId, (*Coeffs_dip)[280][1][1]);
//		printf("\n ProcId %d -- (*Coeffs_dip)[281][1][1]: %f \n", procId, (*Coeffs_dip)[281][1][1]);
//		printf("\n ProcId %d -- (*Coeffs_dip)[560][1][1]: %f \n", procId, (*Coeffs_dip)[560][1][1]);
	}

	return(0);
}


int okadaCoeff2DCFS(float ***Coeffs_st, float ***Coeffs_d, struct pscmp DCFS, struct eqkfm *eqkfm1, struct crust crst, double *strikeR, double *dipR, int full_tensor){

	// [Fahad] Variables used for MPI.
	int procId = 0, numProcs = 1;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	#endif

	double strike, dip, rake;
	double alpha;
	int p1, p, j;
	int Nsel;
	int NF=DCFS.NF;
	int err=0, errp=0;

	// todo [coverage] this block is never tested
	if ((DCFS.nsel!=(*eqkfm1).nsel)){
		print_screen("**Warning: DCFS.nsel!=eqkfm.nsel in okadaCoeff2DCFS. Will use those from eqkfm1. **\n");
		print_logfile("**Error: DCFS.nsel!=eqkfm1.nsel (%d, %d) in okadaCoeff2DCFS.**\n", DCFS.nsel, (*eqkfm1).nsel);
		DCFS.nsel=(*eqkfm1).nsel;
		DCFS.which_pts=(*eqkfm1).selpoints;
	}

	Nsel=DCFS.nsel;

	alpha = (crst.lambda + crst.mu)/(crst.lambda + 2*crst.mu);

	for (int i0=1; i0<=DCFS.nsel; i0++) {
		DCFS.cmb[i0]=0.0;
		for (int a=1; a<=3; a++){
				for (int b=1; b<=3; b++) DCFS.S[i0][a][b]=0;
		}
	}

	#pragma omp parallel for private(p1, p, j, strike, dip, rake) reduction(+:errp)
	for (int i=1; i<=Nsel; i++){

		//-----------------------------------------------------------------------------------------//
		//-----------Calculate Coulomb stress vector for each patch and add them up.---------------//
		//-----------------------------------------------------------------------------------------//

		p1=0;
		for (j=0; j<NF; j++){

			if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
				print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
				print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
				errp+=1;
			}

			for (p=1; p<=eqkfm1[j].np_di*eqkfm1[j].np_st; p++){
				p1+=1;
				DCFS.S[i][1][1]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][1]+eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][1];
				DCFS.S[i][1][2]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][4]+eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][4];
				DCFS.S[i][1][3]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][6]+eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][6];
				DCFS.S[i][2][2]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][2]+eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][2];
				DCFS.S[i][2][3]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][5]+eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][5];
				DCFS.S[i][3][3]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][3]+eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][3];
			}
		}
		DCFS.S[i][2][1]=DCFS.S[i][1][2];
		DCFS.S[i][3][2]=DCFS.S[i][2][3];
		DCFS.S[i][3][1]=DCFS.S[i][1][3];

	}

	//-----------------------------------------------------------------------------------------//
	//-----------------Fill in Syx, Szy, Szx and resolve stress on required plane.-------------//
	//-----------------------------------------------------------------------------------------//

	if (!full_tensor) resolve_DCFS(DCFS, crst, strikeR, dipR, NULL, 1);	//todo not hardwire optrake?

	return(errp!=0);
}

// todo [coverage] this block is never tested
int okadaCoeff_resolve(struct Coeff_LinkList Coeffs, float ***Coeffs_st2, float ***Coeffs_di2, struct crust crst, double *strikeRs, double *dipRs, double *rakeRs){
	//NB: does not have option optrake.

	double *sigma0s;
	double **stress0s;
	double **n, **s, **S, ***S3;
	double strikeR, dipR, rakeR;
	int NTt=omp_get_max_threads(), NT, fm;
	int no_fm_zones=crst.nofmzones;
	int *fm_zones_indices=crst.fmzone;

	n=malloc((no_fm_zones)*sizeof(double *));
	s=malloc((no_fm_zones)*sizeof(double *));
	stress0s=malloc((no_fm_zones)*sizeof(double *));
	sigma0s=dvector(0,no_fm_zones-1);
    S3=d3tensor(0,NTt,1,3,1,3);

	for (int i=0; i<no_fm_zones; i++){

		strikeR=strikeRs[i];
		dipR=dipRs[i];
		rakeR=rakeRs[i];

		n[i]=normal_vector(strikeR, dipR);
		s[i]=slip_vector(strikeR, dipR, rakeR);

		//calculated background stress and normal stress:
		stress0s[i]=mtimesv(crst.S,n[i],NULL,3,3);
		sigma0s[i]=-1.0*vdotv(stress0s[i],n[i],3);		// compression is positive (rock mechanics convention).
	}

	*Coeffs_st2=matrix(1,Coeffs.NP, 1,Coeffs.NgridT);
	*Coeffs_di2=matrix(1,Coeffs.NP, 1,Coeffs.NgridT);

	#pragma omp parallel for private(S, NT, fm)
	for (int i=1; i<=Coeffs.NgridT; i++){
		fm = (no_fm_zones==1) ? 0 :  fm_zones_indices[i];
		NT=omp_get_thread_num();
		S=S3[NT];		//todo test this;

		for (int p1=1; p1<= Coeffs.NP; p1++){

			// find strike coefficients:
			S[1][1]=(double) Coeffs.Coeffs_st[p1][i][1];
			S[1][2]=(double) Coeffs.Coeffs_st[p1][i][4];
			S[1][3]=(double) Coeffs.Coeffs_st[p1][i][6];
			S[2][1]=(double) Coeffs.Coeffs_st[p1][i][4];
			S[2][2]=(double) Coeffs.Coeffs_st[p1][i][2];
			S[2][3]=(double) Coeffs.Coeffs_st[p1][i][5];
			S[3][1]=(double) Coeffs.Coeffs_st[p1][i][6];
			S[3][2]=(double) Coeffs.Coeffs_st[p1][i][5];
			S[3][3]=(double) Coeffs.Coeffs_st[p1][i][3];

			(*Coeffs_st2)[p1][i]= (float) resolve_n(S, n[fm], NULL, crst.fric, stress0s[fm], sigma0s[fm], s[fm]);

			// find dip coefficients:
			S[1][1]=(double) Coeffs.Coeffs_dip[p1][i][1];
			S[1][2]=(double) Coeffs.Coeffs_dip[p1][i][4];
			S[1][3]=(double) Coeffs.Coeffs_dip[p1][i][6];
			S[2][1]=(double) Coeffs.Coeffs_dip[p1][i][4];
			S[2][2]=(double) Coeffs.Coeffs_dip[p1][i][2];
			S[2][3]=(double) Coeffs.Coeffs_dip[p1][i][5];
			S[3][1]=(double) Coeffs.Coeffs_dip[p1][i][6];
			S[3][2]=(double) Coeffs.Coeffs_dip[p1][i][5];
			S[3][3]=(double) Coeffs.Coeffs_dip[p1][i][3];

			(*Coeffs_di2)[p1][i]= (float) resolve_n(S, n[fm], NULL, crst.fric, stress0s[fm], sigma0s[fm], s[fm]);

		}
	}

	free_d3tensor(S3,0,NTt,1,3,1,3);
	for (int i=0; i<no_fm_zones; i++){
		if (s[i]) free_dvector(s[i],1,3);
		free_dvector(n[i],1,3);
		free_dvector(stress0s[i],1,3);
	}
	free(s);
	free(n);
	free(stress0s);
	free_dvector(sigma0s,0,1);
	return(0);
}

// todo [coverage] this block is never tested
int resolvedCoeff2DCFS(float **Coeffs_st, float **Coeffs_d, struct pscmp DCFS, struct eqkfm *eqkfm1, struct crust crst){

	double MaxDCFS=DCFS_cap;
	int p1, p, j;
	int Nsel;
	int NF=DCFS.NF;

	if ((DCFS.nsel!=(*eqkfm1).nsel)){
		print_logfile("**Error: DCFS.nsel!=eqkfm.nsel in resolvedCoeff2DCFS.**\n");
		print_screen("**Warning: DCFS.nsel!=eqkfm.nsel in resolvedCoeff2DCFS. Will use those from eqkfm1. **\n");
		DCFS.nsel=(*eqkfm1).nsel;
		DCFS.which_pts=(*eqkfm1).selpoints;
	}

	Nsel=DCFS.nsel;

	#pragma omp parallel for private(p1, p, j)
	for (int i=1; i<=Nsel; i++){

		//-----------------------------------------------------------------------------------------//
		//-----------Calculate Coulomb stress vector for each patch and add them up.---------------//
		//-----------------------------------------------------------------------------------------//

		//printf("Calculating Okada solutions...");
		DCFS.cmb[i]=0.0;
		p1=0;
		for (j=0; j<NF; j++){
			for (p=1; p<=eqkfm1[j].np_di*eqkfm1[j].np_st; p++){
				p1+=1;
				DCFS.cmb[i]+=eqkfm1[j].slip_str[p]*((double) Coeffs_st[p1][i])+eqkfm1[j].slip_dip[p]*((double) Coeffs_d[p1][i]);
			}
		}
		if (DCFS.cmb[i]>MaxDCFS) DCFS.cmb[i]=MaxDCFS;
		if (DCFS.cmb[i]<-MaxDCFS) DCFS.cmb[i]=-MaxDCFS;
	}
	return 0;
}

int isoDCFS(struct pscmp DCFS, struct eqkfm eqkfm1){
	double M0, r;
	int Nsel;
	double DCFSmax=DCFS_cap;

	// todo [coverage] this block is never tested
	if (DCFS.nsel!=eqkfm1.nsel){
		print_logfile("**Error: DCFS.nsel!=eqkfm.nsel in isoDCFS.**\n");
		print_screen("**Warning: DCFS.nsel!=eqkfm.nsel in isoDCFS. Will choose the one from DCFS.**\n");
		eqkfm1.nsel=DCFS.nsel;
		eqkfm1.selpoints=DCFS.which_pts;
	}

	Nsel=DCFS.nsel;
	#pragma omp parallel for private(r, M0)
	for (int i=1; i<=Nsel; i++){
		r=DCFS.fdist[i];	//r is in km.
		M0=pow(10,1.5*(eqkfm1.mag+6.0));
		DCFS.cmb[i]=(M0/(6.0*PI))*pow(1000*r,-3.0);
		if (fabs(DCFS.cmb[i])>DCFSmax) DCFS.cmb[i]=(DCFS.cmb[i]>0)? DCFSmax : -DCFSmax;
	}
	return(0);
}

//---------------------------------------------------------------------//
//-----					Auxiliary functions						  -----//
//---------------------------------------------------------------------//


int choose_focmec(struct eqkfm eqkfm1, double *strike, double *dip, double *rake){
//assign strike, dip, rake correct value from eqkfm, and convert to radians.

	switch (eqkfm1.whichfm){
		case 1:
			if (strike) *strike=DEG2RAD*eqkfm1.str1;
			if (dip) *dip=DEG2RAD*eqkfm1.dip1;
			if (rake) *rake=DEG2RAD*eqkfm1.rake1;
			break;
		case 2:
			if (strike) *strike=DEG2RAD*eqkfm1.str2;
			if (dip) *dip=DEG2RAD*eqkfm1.dip2;
			if (rake) *rake=DEG2RAD*eqkfm1.rake2;
			break;
		case 0:
			if (strike) *strike=DEG2RAD*eqkfm1.str1;
			if (dip) *dip=DEG2RAD*eqkfm1.dip1;
			if (rake) *rake=DEG2RAD*eqkfm1.rake1;
			break;
		default:
			return(1);
	}
	return (0);
}

void patch_pos(struct eqkfm eqfm, int p, double *east, double *north, double *depth){
//find position of patch in local cartesians.

	double pos_top, dz, dip, strike;

	choose_focmec(eqfm, &strike, &dip, NULL);

	pos_top=eqfm.pos_d[p]-0.5*eqfm.W*(1.0/eqfm.np_di);	//position of top of the patch (along dip).
	*east=eqfm.x+eqfm.pos_s[p]*sin(strike)+pos_top*cos(dip)*cos(strike);
	*north=eqfm.y+eqfm.pos_s[p]*cos(strike)-pos_top*cos(dip)*sin(strike);
	dz=sin(dip)*pos_top;	//depth of top of the patch.
	*depth=eqfm.depth+dz;

	return;
}

// todo [coverage] this block is never tested
double ** comp2tensor(float *v, double ***S0){
	//fills in S0 if given, if NULL allocate new memory.
	double **S;

	S= (S0)? *S0 : dmatrix(1,3,1,3);
	S[1][1]=(double) v[1];
	S[1][2]=(double) v[4];
	S[1][3]=(double) v[6];
	S[2][1]=(double) v[4];
	S[2][2]=(double) v[2];
	S[2][3]=(double) v[5];
	S[3][1]=(double) v[6];
	S[3][2]=(double) v[5];
	S[3][3]=(double) v[3];

	return S;
}

double *normal_vector(double strikeR, double dipR){

	double *n;

	n=dvector(1,3);
	n[1]=-sin(strikeR*DEG2RAD)*sin(dipR*DEG2RAD);
	n[2]=cos(strikeR*DEG2RAD)*sin(dipR*DEG2RAD);
	n[3]=-cos(dipR*DEG2RAD);

	return n;
}

// todo [coverage] this block is never tested
double *slip_vector(double strikeR, double dipR, double rakeR){

	double *s;
	s=dvector(1,3);
	s[1]=cos(rakeR*DEG2RAD)*cos(strikeR*DEG2RAD)+sin(rakeR*DEG2RAD)*sin(strikeR*DEG2RAD)*cos(dipR*DEG2RAD);
	s[2]=cos(rakeR*DEG2RAD)*sin(strikeR*DEG2RAD)-sin(rakeR*DEG2RAD)*cos(strikeR*DEG2RAD)*cos(dipR*DEG2RAD);
	s[3]=-sin(rakeR*DEG2RAD)*sin(dipR*DEG2RAD);

	return s;
}

double *opt_s(double *stress, double sigma, double *n, double *result){
	/* returns direction of maximum shear stress on a plane normal to vector n, subject to stress "stress" and normal stress sigma.
	 * if result ==NULL, ignore it and allocate new memory. Otherwise, memory should have been allocated previously.
	 *
	 * stress[1,2,3]= stress vector;
	 * sigma=normal stress (compression is positive);
	 * n=normal vector to plane;
	 * */

	double s_temp[4], *s;
	double s_tempM;

	s= (result) ? result : dvector(1,3);

	for (int g=1; g<4; g++) s_temp[g]=stress[g]+sigma*n[g];
	s_tempM=norm(s_temp,3);
	for (int g=1; g<4; g++) s[g]= (s_tempM==0)? s_temp[g] : s_temp[g]*(1.0/s_tempM);

	return s;
}

double *sum_v(double *v1, double *v2, double *sum, int N){
// sums two vectors [1...N].
//if sum==NULL, it is ignored and new memory is allocated.

	double *v3;

	if (!v2) return v1;
	if (!v1) return v2;

	v3= (sum) ? sum : dvector(1,N);
	if (v1 && v2) {
		for (int n=1; n<=N; n++) v3[n]=v1[n]+v2[n];
	}
	else {
		if (!v1) v3=v2;
		if (!v2) v3=v1;
	}
	return v3;
}

// todo [coverage] this block is never tested
double resolve_S(double **S, double strikeR, double dipR, double rakeR, double f, double *stress0, double sigma0, double *newrake, int opt_rake){
//To be called after filling in DCFS.S. (with functions below).
//total stress should be passed if optimal rake is to be used (opt_rake==1).
//resolves stress tensor S on mechanisms with dipR, strikeR and rake.
//sigma0, stress0 are not calculated here for efficiency reasons (they don't change as often as S).

	double *n, *s;
	double cmb;

	n=normal_vector(strikeR, dipR);
	if (opt_rake) s=NULL;
	else s=slip_vector(strikeR, dipR, rakeR);

	cmb=resolve_n(S, n, newrake, f, stress0, sigma0, s);

	free_dvector(n,1,3);
	if (s) free_dvector(s,1,3);

	return cmb;
}

double resolve_n(double **S, double *n, double *rake, double fric, double *stress0, double sigma0, double *slip_v){
//resolves stress S on focal mechanism with plane normal to vector n[1..3], and slip direction vector s[1..3].
//if slip_v==NULL, will use optimal rake. This is calculated from stress tensor S as well as background stresses stress0, and background normal pressure sigma0.
//is sigma0==0, stress0==NULL, no background stress.
//sigma0, stress0 are not calculated here for efficiency reasons (they don't change as often as S).

	double *s, *stress, *stress_tot=NULL;
	double sigma, tau;
	int optrake= !slip_v;

	s=dvector(1,3);
	stress=dvector(1,3);

	mtimesv(S,n,stress,3,3);
	sigma=-1.0*vdotv(stress,n,3);	// compression is positive (rock mechanics convention).

	if (optrake) {
		opt_s(stress_tot=sum_v(stress, stress0, NULL, 3), sigma+sigma0, n, s);
		tau=vdotv(stress,s,3);
	}
	else tau=vdotv(stress,slip_v,3);

	if (optrake && rake) (*rake)=RAD2DEG*asin(-s[3]/sqrt(1-n[3]*n[3]));	//todo check if this is ambiguous...

	free_dvector(s,1,3);
	free_dvector(stress,1,3);
	if (optrake) free_dvector(stress_tot,1,3);

	return tau-fric*sigma;

}


