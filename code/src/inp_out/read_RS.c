/*
 * readcatalog.c
 *
 *  Created on: Dec 22, 2011
 *      Author: camcat
 */

#include "read_RS.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int gridPMax=1000;	// max no. points associated with event.
