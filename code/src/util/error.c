/*
 * error.c
 *
 *  Created on: Sep 24, 2013
 *      Author: camcat
 */

#include "error.h"

void error_quit(char * message){
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	if(procId == 0) {
		if (flog){
			fprintf(flog, message);
			fflush(flog);
		}
		printf("%s\n", message);
		fflush(stdout);
	}

	exit(EXIT_FAILURE);
}
