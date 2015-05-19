/*
 * merge.h
 *
 *  Created on: May 18, 2015
 *      Author: camcat
 */


void merge(double a[], int m, double b[], int n, double **sorted, int **ai, int **bi);
void merge_multiple(double **vs, int *lens, int N, double **sorted, int *len_fin, int ***indices);
