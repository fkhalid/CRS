/*
 * Synthetic_stresstests.c
 *
 *  Created on: Apr 17, 2014
 *      Author: camcat
 */

#include <math.h>

#include "defines.h"
#include "inp_out/read_eqkfm.h"


///--------------------------TODO make sure fault is not treated as blind-------------------------//



void Synthetic_stresstests(){

	struct crust crst; crst.mu=26624; crst.lambda=31226; crst.fric=0.3;
	double alpha=(crst.lambda + crst.mu)/(crst.lambda + 2*crst.mu);
	double Sxx, Syy, Szz, Sxy, Syz, Sxz;
	double **S=dmatrix(1,3,1,3);
	double stress;
	double *n, *s;
	struct eqkfm eqkfmco, eqkfmpo;
	double res=5.0;
	int taper=1;	//flag.

	int ndips=3;

	//double dip_slab[21]={10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
	double dip_slab[3]={10, 15, 20};
	double strike_slab=0.0;
	int taper_co[]={0, 1, 1, 1, 1};
	int taper_po[]={0, 1, 1, 1, 1};

	int nWratios=3;
	double Wratios [3]={0.2, 0.5, 1.0};
	double co_depthext=40.0;
	double centroid_depth_co=0.5*co_depthext;

	int nrecfaults=18;
	int nrecfaultdep=3;
	double recfaultdep[3]={0, 5, 10};
	double recfaults_str[18]={-10, 0, 10, -10, 0, 10,-10, 0, 10, 170, 180, 190, 170, 180, 190, 170, 180, 190};
	double recfaults_dip[18]={60, 60, 60, 75, 75, 75, 90, 90, 90, 60, 60, 60, 75, 75, 75, 90, 90, 90};

	int nrec_fault_loc=5;
	//double rec_fault_loc[9]={-50, -40, -30, -20, 0, 20, 30, 40, 50};
	double nrec_fault_area, dx;


	double dmagco=0.05;
	double dmagpo=0.1;
	double magco_min=8.1;
	double magco_max=9.0;
	double magpo_min=6.5;
	double magpo_max=9.0;
	int nmagco=1+(magco_max-magco_min)/dmagco;
	int nmagpo=1+(magpo_max-magpo_min)/dmagpo;
	double *magsco=dvector(0,nmagco-1);
	double *magspo=dvector(0,nmagpo-1);
	for (int m=0; m<nmagco; m++) magsco[m]=magco_min+m*dmagco;
	for (int m=0; m<nmagpo; m++) magspo[m]=magpo_min+m*dmagpo;

	double Wco, Wpo, WR, centr_dist;
	int cc=0, cc0=0;

	FILE *fout0=fopen("summary3.dat","w"), *fout1=fopen("co_stress3.dat", "w"), *fout2=fopen("po_stress3.dat", "w");

	eqkfmco.str1=strike_slab;
	eqkfmco.lat=0.0;
	eqkfmco.lon=0.0;
	eqkfmco.rake1=90;
	eqkfmco.whichfm=1;
	eqkfmco.cuts_surf=1;

	eqkfmco.lat=0.0;
	eqkfmpo.str1=strike_slab;
	eqkfmpo.rake1=90;
	eqkfmpo.whichfm=1;
	eqkfmco.cuts_surf=0;

	eqkfmco.taper=taper_co;
	eqkfmpo.taper=taper_po;


	for (int i=0; i<nmagco; i++){
//	for (int i=nmag-1; i>=0; i--){

		eqkfmco.mag=magsco[i];
		for (int ii=0; ii<nmagpo; ii++){

			eqkfmpo.mag=magspo[ii];

			for (int dip_n=0; dip_n<ndips; dip_n++){
				for (int Wratio=0; Wratio<nWratios; Wratio++){

					cc0++;

					WR=Wratios[Wratio];
					eqkfmco.dip1=dip_slab[dip_n];
					eqkfmco.depth=centroid_depth_co+100;	//+100 to avoid focmec to shift it downwards.
					eqkfmpo.depth=centroid_depth_co+100;	//+100 to avoid focmec to shift it downwards.

					eqkfmpo.dip1=dip_slab[dip_n];
					//eqkfmpo.depth=centroid_depth_po+100;

					Wco=co_depthext/sin(PI*eqkfmco.dip1/180);
					Wpo=Wco*WR;

					focmec2slipmodel(crst, &eqkfmco, res, 0, 0);
					//eqkfmco.L=(eqkfmco.W*eqkfmco.L)/Wco;

					for (int i=1; i<=eqkfmco.np_di*eqkfmco.np_st; i++) {
						eqkfmco.slip_dip[i]=eqkfmco.slip_dip[i]*(eqkfmco.W*eqkfmco.L)/(Wco*500);
						eqkfmco.pos_d[i]=eqkfmco.pos_d[i]*Wco/eqkfmco.W;
					}
					eqkfmco.L=500;
					eqkfmco.W=Wco;
					eqkfmco.depth=centroid_depth_co;
					if (taper){
						suomod1_resample(eqkfmco, &eqkfmco, res, 0.0);	//create a slip model with right resolution.
						suomod1_taper(eqkfmco, &eqkfmco);
					}

					//centr_dist=0.5*(co_depthext+po_depthext)*tan(PI*eqkfmco.dip1/180);
					centr_dist=fabs(0.5*co_depthext*(1+WR))/tan(PI*eqkfmco.dip1/180);
					eqkfmpo.lon=(centr_dist/Re)*(180/pi);
					focmec2slipmodel(crst, &eqkfmpo, res, 0, 0);
					//rescale slip to required area:
					for (int i=1; i<=eqkfmpo.np_di*eqkfmpo.np_st; i++) {
						eqkfmpo.slip_dip[i]=eqkfmpo.slip_dip[i]*(eqkfmpo.W*eqkfmpo.L)/(WR*eqkfmco.W*eqkfmco.L);
						eqkfmpo.pos_d[i]=eqkfmpo.pos_d[i]*WR*eqkfmco.W/eqkfmpo.W;
					}
					eqkfmpo.depth=co_depthext*(1.0+0.5*WR);
					eqkfmpo.L=eqkfmco.L;
					eqkfmpo.W=WR*eqkfmco.W;
					eqkfmco.slip_str[1]=eqkfmpo.slip_str[1]=0.0;
					if (taper){
						suomod1_resample(eqkfmpo, &eqkfmpo, res, 0.0);	//create a slip model with right resolution.
						suomod1_taper(eqkfmpo, &eqkfmpo);
					}
					nrec_fault_area=eqkfmpo.W*cos(PI*eqkfmpo.dip1/180);	//area of faults above afterslip

//					char fname[120];
//					sprintf(fname,"coseismic_models/co%.03d", cc0);
//					print_slipmodel(fname, &eqkfmco, 1);
//
//					sprintf(fname,"postseismic_models/po%.03d", cc0);
//					print_slipmodel(fname, &eqkfmpo, 1);

					printf("Centr_dist=%.5e]\n", centr_dist);
					printf("Co: [mag=%.3lf\tW=%.3lf\tL=%.3lf\tS=%.3e\tNP=%d]\n", eqkfmco.mag, eqkfmco.W, eqkfmco.L, eqkfmco.slip_dip[1], eqkfmco.np_di*eqkfmco.np_st);
					printf("Po: [mag=%.3lf\tW=%.3lf\tL=%.3lf\tS=%.3e\tNP=%d]\n\n", eqkfmpo.mag, eqkfmpo.W, eqkfmpo.L, eqkfmpo.slip_dip[1], eqkfmpo.np_di*eqkfmpo.np_st);
					fflush(stdout);

					for (int rfl=0; rfl<nrec_fault_loc; rfl++){
							for (int rfl2=0; rfl2<nrecfaultdep; rfl2++){

							//dx=(nrec_fault_loc==1) ? 0 : -0.5*nrec_fault_area+rfl*nrec_fault_area/(nrec_fault_loc-1);	//on afterslip
							dx=(nrec_fault_loc==1) ? 0 : -0.5*nrec_fault_area+3*rfl*nrec_fault_area/(nrec_fault_loc-1);	//further inland

							pscokada(0, 0, eqkfmco.depth,  eqkfmco.str1,  eqkfmco.dip1, eqkfmco.L, eqkfmco.W,  0.0, -1.0*eqkfmco.slip_dip[1],
									0.0, centr_dist+dx, recfaultdep[rfl2], &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz, alpha, crst.lambda, crst.mu, crst.fric);
									S[1][1]=1e6*Sxx;
									S[2][2]=1e6*Syy;
									S[3][3]=1e6*Szz;
									S[1][2]=1e6*Sxy;
									S[2][3]=1e6*Syz;
									S[1][3]=1e6*Sxz;


							for (int rf=0; rf<nrecfaults; rf++){
								n=normal_vector(recfaults_str[rf], recfaults_dip[rf]);
								s=slip_vector(recfaults_str[rf], recfaults_dip[rf], 0.0);
								stress=resolve_n(S, n, NULL, crst.fric, NULL, 0.0, s);
								fprintf(fout1, "%.5e\t", stress);
							}

							pscokada(0, centr_dist, eqkfmpo.depth,  eqkfmpo.str1,  eqkfmpo.dip1, eqkfmpo.L, eqkfmpo.W,  0.0, -1.0*eqkfmpo.slip_dip[1],
									0.0, centr_dist+dx, recfaultdep[rfl2], &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz, alpha, crst.lambda, crst.mu, crst.fric);
									S[1][1]=1e6*Sxx;
									S[2][2]=1e6*Syy;
									S[3][3]=1e6*Szz;
									S[1][2]=1e6*Sxy;
									S[2][3]=1e6*Syz;
									S[1][3]=1e6*Sxz;


							for (int rf=0; rf<nrecfaults; rf++){
								n=normal_vector(recfaults_str[rf], recfaults_dip[rf]);
								s=slip_vector(recfaults_str[rf], recfaults_dip[rf], 0.0);
								stress=resolve_n(S, n, NULL, crst.fric, NULL, 0.0, s);
								fprintf(fout2, "%.5e\t", stress);
							}

							for (int rf=0; rf<nrecfaults; rf++){
								cc++;
								fprintf(fout0, "%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.3f\n", eqkfmco.mag, eqkfmpo.mag, eqkfmco.slip_dip[1], eqkfmpo.slip_dip[1], dip_slab[dip_n], Wratios[Wratio], recfaults_str[rf], recfaults_dip[rf], centr_dist, dx/nrec_fault_area);
							}
						}
					}
					//printf("\n");
				}
			}
			fprintf(fout2,"\n");
			fprintf(fout1,"\n");
		}
		//printf("\n-------------------------------------------------\n");
		fprintf(fout1,"\n");
		fprintf(fout2,"\n");

	}
}
