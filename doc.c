/*
 * doc.c
 *
 *  Created on: Jun 25, 2013
 *      Author: camcat
 */

//------------------------------------------------------------------------------------------|
//	NAME	|	MEANING																		|
//------------------------------------------------------------------------------------------|
	Nm		|	Number of mainshocks (events with slip models available).
	Nas		|	Number of afterslip snapshots.
	NFtot	|	Total numbner of faults of mainshocks (sum of Nfaults).
	Ntot	|	Number of aftershocks used as sources.
	Nfm		|	Number of sources with foc. mec.
	L		|	Number of time steps fr integration (and at which afterslip is calculated).
	NgridT	|	grid points.


//----------------------------------------------------------------------------------------------------------------------|
//	TYPE	|	VARIABLE		|		INDICES			|	CONTENT														|
//----------------------------------------------------------------------------------------------------------------------|
	char **	|	slipmodels		|		0...Nm-1		|	file names of slip models for mainshocks.
	char ** | 	multimodels		|		0...Nm-1		|	file names of multiple slip models for mainshocks.
	int *   | 	no_slipmodels	|		0...Nm-1		|	number of slip models available for each mainshock.
	double*	|	tmain			|		0...Nm-1		|	time of mainshocks.
	int *	|	mmain			|		0...Nm-1		|	magnitude of mainshocks.
	int*	|	Nfaults			|		0...Nm-1		|	number of faults for each mainshock.
	int *	|	which_main		|		0...Nm-1		|	elements of eqkfm1 to which mainshocks correspond.

	char **	|	afterslipmodels	|		0...Nas-1		|	file names of afterslip snapshots.
	double*	|	t_afterslip		|		0...Nas-1		|	times of afterslip snapshots.
	eqkfm * | 	eqkfm_aft0		|	0...Nfaults[0]-1	|	temporary structure into which afterslip files are red.
	eqkfm *	| 	eqkfm_aft		|  0...Nas*Nfaults[0]-1	|	resampled afterslip snapshots.
	eqkfm *	|	eqkfm_aftsplines|	0...Nfaults[0]*(L+1)|	interpolated afterslip snapshots (non cumulative). Last Nfaults[0] elements: cumulative slip.
	eqkfm *	|	eqkfm0			|		0...NFtot-1		|
	eqkfm *	|	eqkfm0res		|		0...NFtot-1		|
	eqkfm *	|	eqkfm1			|		0...Ntot-1		|	Aftershocks used as sources.
	eqkfm *	|	eqkfm1fm		|		0...Nfm-1		|	Aftershocks used as sources, from foc. mec. catalog.
	pscmp *	|	DCFS			|  0...Ntot-1 (Nm-1)	|	discrete stress changes with (without) including aftershocks.

	double*	|	times2			|		0...L-1	(or L?)	|	time steps for continuous process (afterslip).
	double*	|	tevol_afterslip	|		0...L-1			|	afterslip coeffiecients at each time step (no splines).
	doub**	|	DCFSrand		|	0...L-1; 1...NgridT	|	perturbed afterslip field.

//	----------------------------------//
//	in forecast_stepG2_new.c:
//----------------------------------------------------------------------------------------------------------------------|
//	TYPE	|	VARIABLE		|		INDICES			|	CONTENT														|
//----------------------------------------------------------------------------------------------------------------------|

	sh.int *|	num_eqk			|		1...N			|	number of earthquakes affecting each grid point.
	sh.int *|	which_eqk		|	1...N; 0...Neq		|	list of eqks affecting each grid point.
	sh.int *|	which_elem		|	1...N; 0...Neq		|	index of gridpoint in structure DCFS of eqks affecting it.

