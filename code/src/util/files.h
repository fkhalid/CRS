/******************************************************************************
  "files.c"

  description:  Helper for working with files


  project:  Simulation and parameter estimation using ETAS and shakemaps
   author:
     date:
******************************************************************************/


int fgetword(FILE *fp, char *word, int count);
/* Liest ein von Trennzeichen begrenztes Wort ein */

char *sgetint(char *, int *);
/* liest aus String s eine von Blanks begrenzte int-Zahl ein */

int fgetsequ(FILE *, char *, int);
/*  liest eine Sequenz bestimmter Laenge ein */

int fgetline(FILE *, char *, int);
/* liest eine ganze Zeile ein */

int isalnum_0(int);
/* prueft, ob Zeichen alphanumerisch. Null ist Begrenzer! */

int isalnumde(int);
/* prueft, ob Zeichen alphanumerisch im erweiterten Zeichensatz */
