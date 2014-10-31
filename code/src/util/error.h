/*
 * error.h
 *
 *  Created on: Sep 24, 2013
 *      Author: camcat
 */

#ifndef ERROR_H_
#define ERROR_H_

void error_quit_fun(char *fun, const char * format, ...);
void print_logfile_fun(char *fun, const char * format, ...);
void print_screen_fun(char *fun, const char * format, ...);

#endif /* ERROR_H_ */
