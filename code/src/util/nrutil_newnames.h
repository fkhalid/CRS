//TODO add header (also to .c file) and change name.

#define NR_END 1
#define FREE_ARG char*

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *);
void nrerror(char error_text[]);
int *iarray(long nl, long nh);
double *darray(long nl, long nh);
unsigned long *larray(long nl, long nh);
double **d2array(long nrl, long nrh, long ncl, long nch);
float **f2array(long nrl, long nrh, long ncl, long nch);
int **i2array(long nrl, long nrh, long ncl, long nch);
float ***f3array(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_iarray(int *v, long nl, long nh);
void free_darray(double *v, long nl, long nh);
void free_larray(unsigned long *v, long nl, long nh);
void free_d2array(double **m, long nrl, long nrh, long ncl, long nch);
void free_f2array(float **m, long nrl, long nrh, long ncl, long nch);
void free_i2array(int **m, long nrl, long nrh, long ncl, long nch);
void free_f3array(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
