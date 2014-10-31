/*
 * propagate_results.c
 *
 *  Created on: Aug 30, 2013
 *      Author: camcat
 */

#include "propagate_results.h"

//int check_if_snapshot_exists(char * folder, double *t,  struct tm reftime, long *hash){
////returns: 0 if file not found, 1 if found, -1 if error.
////if file is found, the time value is read into t;
////if t or hash is null, file will not be read.
//
//	/*simpler version:
//
//	char fname[120];
//	sprintf(fname,"%s/%s",folder,check_if_snapshot_filename);
//	if (fopen(fname, "r") == NULL) return 0;
//	else return 1;
//	 */
//
//	DIR *pdir = NULL;
//	struct tm oldtime;
//	struct dirent *pent = NULL;
//	int Nchar=120;
//	char fname[120], fname2[120], line[120];
//
//	FILE *fin;
//
//	sprintf(fname, check_if_snapshot_filename);
//
//	pdir = opendir (folder); // "." will refer to the current directory
//	if (pdir == NULL)
//	{
//	    printf ("\nERROR! pdir could not be initialised correctly");
//	    return -1;
//	}
//
//	while ((pent = readdir(pdir)))
//	{
//	    if (pent == NULL)
//	    {
//	    	printf ("ERROR! pent could not be initialised correctly");
//	        return -1;
//	    }
//	    if (strcmp (pent->d_name, fname)==0) {
//	    	sprintf(fname2,"%s/%s",folder, fname);
//	    	if (t && hash){
//				fin=fopen(fname2,"r");
//
//				fgets(line,Nchar,fin);
//				sscanf(line, "%d-%d-%dT%d:%d:%dZ", &(oldtime.tm_year), &(oldtime.tm_mon), &(oldtime.tm_mday), &(oldtime.tm_hour), &(oldtime.tm_min), &(oldtime.tm_sec));
//				oldtime.tm_year-=1900;
//				oldtime.tm_mon-=1;
//				oldtime.tm_isdst=0;
//				if (t) *t=difftime(mktime(&oldtime), mktime(&reftime))*SEC2DAY;
//				fgets(line,Nchar,fin);
//				if (hash) sscanf(line, "%ld", hash);
//				fclose(fin);
//	    	}
//	    	return 1;
//	    }
//	}
//
//    closedir (pdir);
//
//    return 0; // file not found.
//}
//
//int load_oldLL(char * folder, double ***LLs){
//	// returns 1 if errors, 0 otherwise.
//
//	int N1, N2=7, NH=1;	//rows, columns, headers.
//	char fname[120];
//
//	sprintf(fname,"%s/%s",folder,LLsnapshot_filename);
//
//	N1=countline(fname);
//	if (N1<=0) return 1;
//	*LLs=dmatrix(1,N2,1,N1-NH);
//
//	read_matrix(fname, N2, NH, *LLs, 0);
//
//	return 0;
//
//}
//
//void load_gammas(char * folder, int p, double **gammas_old, int NG){
//
//	char fname[120];
//	sprintf(fname,"%s/%d_%s",folder,p,gammas_filename);
//
//	int NH=0;	//rows, columns, headers.
//
//	read_matrix_transpose(fname, NG, NH, gammas_old, 0);
//
//}
//
//void write_gammas(char *folder, int p, double **gammas, int Nsur, int NG){
//
//	char fname[120];
//	FILE *fout;
//
//	sprintf(fname,"%s/%d_%s",folder,p,gammas_filename);
//	fout= fopen(fname, "w");
//	for (int ns=1; ns<=Nsur; ns++){
//		for (int n=1; n<=NG; n++) fprintf(fout,"%.3e\t",gammas[ns][n]);
//		fprintf(fout,"\n");
//	}
//	fclose(fout);
//
//}
