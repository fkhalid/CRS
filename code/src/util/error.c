/*
 * error.c
 *
 *  Created on: Sep 24, 2013
 *      Author: camcat
 */

#include "error.h"

void error_quit(char * message){
	if (flog){
		fprintf(flog, message);
		fflush(flog);
	}
	printf("%s\n", message);
	fflush(stdout);
	exit(EXIT_FAILURE);
}
