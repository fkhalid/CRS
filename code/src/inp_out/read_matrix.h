/******************************************************************************
  "read_matric.h"

  description:  Einlesen von Werten aus einer Datei in eine Matrix


  project:  Simulation and parameter estimation using ETAS and shakemaps
   author:  Christoph Bach
     date:  2010-07-27
******************************************************************************/ 
// vormals einlese.h
#ifndef READ_MATRIX_H
#define READ_MATRIX_H

#include <stdio.h>
//#include <string.h>

#include "../defines.h"
#include "../util/files.h"

#define countcol(f) countcol_header(f, 0)

//-----------------------------------------------------------------------------
// INPUT:       infile     :  Name der einzulesenden Datei
//              columns    :  Anzahl der Spalten in dieser Datei
//              headerlines:    "     "  Titelzeilen
// OUTPUT:      data       :  Matrix (S,Z) mit S  Spalten Z  Zeilen
//              rows       :  Zeilenzahl
//-----------------------------------------------------------------------------
int read_matrix(char *infile,int columns, int headerlines, double **data, long *rows);
int read_matrix_transpose(char *infile, int columns, int headerlines, double **data, long *rows);
int read_matrixT(char *infile,int columns, int headerlines, double **data, char**files, long *rows);
int countline(char *filename);
int countcol_header(char *filename, int headerlines);

#endif // READ_MATRIX_H
