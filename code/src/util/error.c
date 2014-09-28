/*
 * error.c
 *
 *  Created on: Sep 24, 2013
 *      Author: camcat
 */

#include "error.h"
#include "../defines.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>


#ifdef _CRS_MPI
	#include "mpi.h"
#endif

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


void print_logfile(const char * format, ...)
{
  char buffer[3000];
  int procId = 0;
  va_list args;
  va_start (args, format);

  	#ifdef _CRS_MPI
  		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  	#endif

  	if(procId == 0) {
  		if (flog){
  			  vsprintf (buffer,format, args);
  			  fprintf(flog,"%s",buffer);
  			  fflush(flog);
  		}
  	}
  	va_end (args);
}

void print_screen(const char * format, ...)
{
  char buffer[3000];
  int procId = 0;
  va_list args;
  va_start (args, format);

  	#ifdef _CRS_MPI
  		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  	#endif

  	if(procId == 0) {
  		if (verbose_level){
		  vsprintf (buffer,format, args);
  		  printf("%s",buffer);
  		  fflush(stdout);
  		}
  	}
  	va_end (args);
}
