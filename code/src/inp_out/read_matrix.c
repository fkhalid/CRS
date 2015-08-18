/******************************************************************************
 "read_matrix.c"

 description:  Einlesen von Werten aus einer Datei in eine Matrix


 project:  Simulation and parameter estimation using ETAS and shakemaps
 author:  Christoph Bach
 date:  2010-07-27
 ******************************************************************************/

#include "read_matrix.h"

#include <stdio.h>

#include "../defines.h"
#include "../util/files.h"

int read_matrix(char *infile, int columns, int headerlines, double **data, long *rows) {
/* Reads an ascii file into a double array.
 *
 * Input:
 *  infile: input file name
 *  columns: no. of columns
 *  headerlines: no. of header lines, to be skipped.
 *
 * Output:
 *  data: contains the values from the file. Memory must be allocated previously.
 *  rows: no. of rows.
 */


	FILE *fin;
	char title[3000];
	double Zd, dum;
	int s, z, i, ans, dumerror;
	long ZZ, N;

	if((fin = fopen(infile, "r"))==NULL){
		print_screen("**Error: unable to open input file %s (inread_matrix).**\n", infile);
		print_logfile("**Error: unable to open input file %s (inread_matrix).**\n", infile);
		return (1);
	}

	for (i = 1; i <= headerlines; i++)
		fgetline(fin, title, 300);

	ans = 1;
	N = 0;

	while (ans != EOF) {
		ans = fscanf(fin, "%lf", &dum);
		if (ans != EOF) N++;
	}

	fclose(fin);
	Zd = N * 1.0 / (1. * columns);
	if (rows) *rows = (long) Zd;	//if null pointer is passed, ignore.
	ZZ = (long) Zd;

	if (Zd - 1.0 * ZZ != 0.0) {
		print_screen("Error: Mismatch in the number of columns!");
		print_logfile("Error: Mismatch in the number of columns!");
		return(1);
	}

	fin = fopen(infile, "r");
	for (i = 1; i <= headerlines; i++)
		fgetline(fin, title, 300);
	for (z = 1; z <= ZZ; z++)
		for (s = 1; s <= columns; s++)
			dumerror = fscanf(fin, "%lf", &data[s][z]);
	fclose(fin);

	return(0);

}

int countline(char *filename){
/* Returns the number of lines in file "filename"
 */

	FILE *fin;
	int dum=0;
	int counter =0;

	if((fin = fopen(filename, "r"))==NULL) {
		print_screen(" **Error: unable to open input file %s. (countline.c)**\n", filename);
		print_logfile(" **Error: unable to open input file %s. (countline.c) **\n", filename);
		return -1;
	}

	while(dum!=EOF) {
	      dum=fgetc(fin);
	      if(dum=='\n')  counter++;
	    }

	fclose(fin);
	return counter;
}

int countcol_header(char *filename, int headerlines){
	/* Returns the number of columns in file "filename", using the first line after headerlines.
	 */

	char title[3000];
	FILE *fin;
	int dum=0;
	int counter =0;

	if((fin = fopen(filename, "r"))==NULL) {
		print_screen(" **Error: unable to open input file %s.**\n", filename);
		print_logfile(" **Error: unable to open input file %s.**\n", filename);
		return -1;
	}

	if (countline(filename)<=headerlines){
		print_screen(" **Error: headerlines should be less than total no. of lines in file %s.**\n", filename);
		print_logfile(" **Error: headerlines should be less than total no. of lines in file %s.**\n", filename);
		return -1;
	}

	for (int i = 1; i <= headerlines; i++){
			fgetline(fin, title, 3000);
	}
	dum=(int) '\t';
	while(dum=='\t' || dum==' ' || dum=='\n') dum=fgetc(fin);

	while (dum!='\n'){
		while(dum!='\t' && dum!=' ' && dum!='\n') {
			dum=fgetc(fin);
		}
		if(dum!='\n')  {
			counter++;
			while(dum=='\t' || dum==' ') dum=fgetc(fin);
			if (dum=='\n') counter-=1;
		}
	}
	fclose(fin);

	return counter+1;
}












