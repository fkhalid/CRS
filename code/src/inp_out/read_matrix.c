/******************************************************************************
 "read_matrix.c"

 description:  Einlesen von Werten aus einer Datei in eine Matrix


 project:  Simulation and parameter estimation using ETAS and shakemaps
 author:  Christoph Bach
 date:  2010-07-27
 ******************************************************************************/
// vormals einlese.c

#include "read_matrix.h"

#include <stdio.h>
//#include <string.h>

#include "../defines.h"
#include "../util/files.h"

//-----------------------------------------------------------------------------
// INPUT:       infile     :  Name der einzulesenden Datei
//              columns    :  Anzahl der Spalten in dieser Datei
//              headerlines:    "     "  Titelzeilen
// OUTPUT:      data       :  Matrix (S,Z) mit S  Spalten Z  Zeilen
//              rows       :  Zeilenzahl
//-----------------------------------------------------------------------------
int read_matrix(char *infile, int columns, int headerlines, double **data, long *rows) {
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

int read_matrix_transpose(char *infile, int columns, int headerlines, double **data, long *rows){

	FILE *fin;
	char title[300];
	double Zd, dum;
	int s, z, i, ans, dumerror;
	long ZZ, N;

	if((fin = fopen(infile, "r"))==NULL){
		print_screen(" **Error: unable to open input file %s. (nread_matrix_transpose)**\n", infile);
		print_logfile(" **Error: unable to open input file %s. (read_matrix_transpose)**\n", infile);
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
			dumerror = fscanf(fin, "%lf", &data[z][s]);
	fclose(fin);

	return(0);


}
int read_matrixT(char *infile, int columns, int headerlines, double **data, char**files, long *rows) {
	FILE *fin;
	char title[300], dum;
	double Zd;
	int s, z, i, ans, dumerror;
	long ZZ, N;

	if((fin = fopen(infile, "r"))==NULL){
		print_screen(" **Error: unable to open input file %s. (nread_matrixT)**\n", infile);
		print_logfile(" **Error: unable to open input file %s. (read_matrixT)**\n", infile);
		return (1);
	}

	for (i = 1; i <= headerlines; i++)
		fgetline(fin, title, 300);

	ans = 1;
	N = 0;


	while (ans != EOF) {
		ans = fscanf(fin, "%s", &dum);
		if (ans != EOF) N++;
	}

	fclose(fin);
//	rewind(fin);

	Zd = N * 1.0 / (1.0 * columns);
	*rows = (long) Zd;
	ZZ = (long) Zd;


	if (Zd - 1.0 * ZZ != 0.0) {
		print_screen("Error: Mismatch in the number of columns!");
		print_logfile("Error: Mismatch in the number of columns!");
		return(1);
	}

	char dums[80];
	fin = fopen(infile, "r");
	for (i = 1; i <= headerlines; i++)
		fgetline(fin, title, 300);
	for (z = 0; z < ZZ; z++){
		for (s = 0; s < columns-1; s++){
			dumerror = fscanf(fin, "%lf", &(data[s+1][z+1]));
		}
		dumerror = fscanf(fin, "%s", files[z]);
	}
	return(0);

}

int countline(char *filename){
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
	/*
	 * Returns the number of columns contained in a file, using the first line after headerlines.
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












