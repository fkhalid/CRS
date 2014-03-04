/*
 * error.c
 *
 *  Created on: Sep 24, 2013
 *      Author: camcat
 */

#include "error.h"

void error_quit(char * message){
	printf("%s\n", message);
	exit(EXIT_FAILURE);
}
