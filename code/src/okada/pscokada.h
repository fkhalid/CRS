

/************************************************************/
/*						cmbfix								*/
/*				calculate Coulomb Stress					*/
/*															*/
/* input:													*/
/*  stress tensor, pore pressure, friction coefficient		*/
/*  rupture orientation parameter (strike, dip and rake)	*/
/*															*/
/* return:													*/
/*  Coulomb stress (cmb) and normal stress (sig)			*/
/************************************************************/

void cmbfix(double sxx, double syy, double szz, double sxy, double syz, double szx, double p, double f, double *cmb, double *sig, double st, double di, double ra);


/************************************************************/
/*						pscokada							*/
/************************************************************/

void pscokada(double x1, double y1, double z1, double strike1, double dip1, double L, double W, double slip_strike, double slip_dip, double open, double x2, double y2, double z2, double *sxx, double *syy, double *szz, double *sxy, double *syz, double *szx, double alpha, double lambda, double mu, double friction);
