/*
 * resolve_onfault.h
 *
 *  Created on: Jun 24, 2014
 *      Author: camcat
 */

#ifndef RESOLVE_ONFAULT_H_
#define RESOLVE_ONFAULT_H_

int read_inputfile2(char *input, char *slipmodelfile, double *lat, double *lon, double *dep, double *str, double *dip, double *rake, int *optrake);
int resolve_onfault(char **argv);

#endif /* RESOLVE_ONFAULT_H_ */
