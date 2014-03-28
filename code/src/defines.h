/******************************************************************************
  "defines.h"

  project:  Simulation and parameter estimation using ETAS and shakemaps
   author:  Christoph Bach
            Sebastian Hainzl
     date:  2010-07-30
******************************************************************************/

#ifndef DEFINES_H
#define DEFINES_H

// ----- [Fahad] Added for MPI -----
#define _CRS_MPI						// FIXME: Should be set depending on whether of not mpicc is used ...
#define BCAST_FLAGS_SIZE 14				// No. of scalar variables in 'struct flags'
#define SIZE_BCAST_MODEL_PARAMETERS 38	// No. of scalar variables in 'struct BCast_Model_Parameters'
// ---------------------------------

#define log10(A) (log(A)/log(10))
#define sq(A) ((A)*(A))
#define sign(A) ((int) (A/fabs(A)))
#define max(A,B) (A>B)? A : B
#define eps0 1.0e-3
#define PI (3.141592653589793)
#define pi (3.141592653589793)
#define EUL (2.718281828459)
#define ASIZE 100000
#define CSIZE 512
#define Re (6370)
#define DEG2RAD  (0.0174532925)
#define KM2M     (1000.0)
#define RAD2DEG  (57.2957795147)
#define SEC2DAY	(1.0/(24.0*3600.0))
#define T2SEC(x) (x/(double)CLOCKS_PER_SEC)
#define tol0 1e-10	//tolerance for double comparison.

#include <stdio.h>

#define check_if_snapshot_filename "snapshot.info"	//file which exists only if at least one forecast has been already been produced, and contains time to which other files gammas and LL) refer.
#define LLsnapshot_filename "LL.dat"	//file which exists only if at least one forecast has been already been produced, and contains time to which other files gammas and LL) refer.
#define gammas_filename "gammas.dat"	//file which exists only if at least one forecast has been already been produced, and contains time to which other files gammas and LL) refer.
#define old_LLfolder "output/oldLL"
#define logfolder "output/log/"

/*	verbose_level: only global variable. settings:
 *
 * 	0: output nothing.
 *  1: output to screen: main operations and errors which lead to abort.
 *  2: output to screen: also minor errors and warnings. (also writes log file).
 *  3: output modified slip models, stress fields etc. (produces several files!)
 *  4: output even more files (e.g. Fourier values etc...)
 *
 * */	//todo implement these!!

extern char cmb_format[120];
extern int verbose_level;
extern int gridPMax;	//todo global variables are bad...
extern double DCFS_cap;
extern FILE *flog;

struct flags{
	int err_recfault;
	int err_slipmodel;
	int err_afterslipmodel;
	int err_gridpoints;
	int OOPs;
	//afterslip:
	int afterslip;
	int splines;
	//control way aftershocks are treated:
	int aftershocks;
	int only_aftershocks_withfm;
	int full_field;			  //if (2): full field for all events.  (1) use available foc mec. (0) use isotropic field for all.
	int aftershocks_fixedmec; //controls is fixed foc. mec. should be used for events w/o foc mec, when fullfield=2 (otherwise, will draw a random one).
	int aftershocks_mode;
	//these can change with each iteration:
	int new_slipmodel;
	int sample_all;
};

struct set_of_models{
	int Nmod;
	int *NF_models;	//no. of faults for each model;
	int NFmax;
	int current_model;
	struct eqkfm *set_of_eqkfm; //contains all models. indices: [0...sum(NF_models)-1].
};

//Linked list with okada coefficients between fault patches and cells. Each element represents one earthquake (also with multiple faults).
struct Coeff_LinkList{
	int NF;				// tot. no of faults;
	int which_main;		// index of pscmp DCFS to which earthquake refer;
	int NP;				// tot. no. of patches (sum of no. of patches of individual faults);
	int NgridT;			// no. of grid cells.
	float ***Coeffs_st, ***Coeffs_dip;	// Coefficient for strike slip, dip slip displacements.
	struct Coeff_LinkList *next;	// pointer to next element.
};

// earthquake catalog.
struct catalog{
	//properties of earthquakes:
	double Mc, b;
	double *t;			//time
	double *mag;		//magnitude;
	double *lat0;
	double *lon0;
	double *x0;
	double *y0;
	double *depths0;
	double *err;
	double *verr;
	int *ngrid;			//no. of of cells associated with each earthquake
	int **ngridpoints;	//indices of cells associated with each earthquake:	ngridpoints[eqk][cell_index].
	double **weights;	//weight of cells associated with each earthquake:	weights[eqk][cell_weight].
	struct crust *pcrst;
	//general catalog properties:
	long Z;
	double tstart;
	double tend;
};

struct pscmp{
	double m;			// event magnitude
	double t;			// event time
	double *fdist;		// distance to fault
	double ***S;		// stress tensor S[cell_index][i][j].	cell_index: 1...nsel.
	double ***S1;		// also stess tensor, may be used if 2 foc. mec. are available.
	double *cmb;		// coulomb stress cmb[cell_index], cell_index: 1...nsel.
	double *cmb0;		// stores undisturbed cmb field (if the only source of errors is the grid point uncertainty).
	double *Dcmb;		// stores range of cmb field (if the only source of errors is the grid point uncertainty).
	long    Z;
	//these may be used to save focal mechanism parameters of OOPs.
	double *st1;
	double *di1;
	double *ra1;
	double *st2;
	double *di2;
	double *ra2;
	//focal mechanism zone index for each grid point:
	int nsel;			//no. of cell points affected by this event.
	int *which_pts;		//indices of cell points affected by this event (relative to arrays lat, lon, depth). range: 1...nsel.
	//int nLat,nLon,nD;	//no. of cells with unique lat, lon, depth (describing overall geometry).
    int index_cat;		//index of event in catalog (only for the catalog used for LL calculation). set to 0 if event is not in catalog.
    int NF; 			//number of faults of mainshock (i.e. no. of eqkfm object mapping to this event).
};

struct crust{
//describes elastic properties of the entire crustal volume.

	//physical properties:
	double str0, dip0, rake0; // orientation of best oriented mechanism in regional stress field: this is redundant given stress tensor below, but easier to keep if for later.
	double fric;		// coefficient of friction.
	double skepton;		// skeption coefficient
	double **S;			// regional stress tensor;
	double lambda;		// lame' parameter
	double mu;			// lame' parameter
	// coords. of domain
	double lat0, lon0;
	double 	latmin, \
			latmax, \
			lonmin, \
			lonmax, \
			depmin, \
			depmax;
	// grid for forecast:
	int nLat_out, \
		nLon_out, \
		nD_out;
	double dlat_out;
	double dlon_out;
	double ddepth_out;
	// grid for calculations:
	int N_allP;			// no. of points (should be same as DCFS0).
	int *list_allP;		// list. of points (should be same as DCFS0): [1,2,3,...N_allP].
	double dlat;		// spacing (lat)
	double dlon;		// ...
	double ddepth;		// ...
	double dmags;
	int nLat;
	int nLon;
	int nD;
	int nmags;
	// coordinates of all grid cell centers.
	double 	*lat, \
			*lon, \
			*depth;
	// coordinates of all grid cell centers (for output).
	double 	*lat_out, \
			*lon_out, \
			*depth_out;
	double *dAgrid;
	double *x;
	double *y;
	double *rate0;	//adds up to 1.
	double r0;		//daily rate for entire region.
	double *mags, *GRmags;	//GRmags=Gutenberg-Richter coefficients corresponding to each magnitude bin.
	int nofmzones;	//no of zones characterized by a different set of receiver faults.
	int *fmzone;		//list of fm zones for each grid point.
	int uniform;
};

struct slipmodels_list{
	int constant_geometry;
	int is_afterslip;
	int NSM;	//no. of events.
	int *Nfaults;
	int *no_slipmodels;
	double *tmain;	//times.
	double *mmain;	//magnitudes.
	double *disc;
	char **slipmodels;
};

//structure describing a single fault earthquake (arrays can be used to describe multiple fault events).
struct eqkfm{	//for events on multiple faults, use a list of these.
	int is_mainshock;	//mainshock with multiple patches (treated differently when calculating DCFS).
	int is_slipmodel;	//if set to 0, focal mechanism is not available.
	int np_st, np_di;	//no. of patches along strile, no. of patches along dip,
	int whichfm;		//index of foc. mec. to use (0=both;1;2).
	int nsel;			//no. of cell points affected by this event.
	int noise;			// flag (if set, tot_slip refers to model from which noise was generated):
	int *taper;		// taper [1...4]=top,bottom,right,left.
	double t;		//time of event.
	double lat;		//0_lat in Wang input file;
	double lon; 	//0_lat in Wang input file;
	double depth; 	//0_lat in Wang input file;
	double x;		//eastwards coordinate in local system;
	double y;		//northward coordinate in local system.
	double mag;		//magnitude
	double tot_slip;//tot. slip on the fault (sum of scalar value of slip for each patch).
	double L;		//fault length
	double W;		//fault width
	//focal mechanism parameters for both possible planes:
	double str1;
	double str2;
	double dip1;
	double dip2;
	double rake1;
	double rake2;
	//vector containing one element per patch:
	//NB: by convention, slip_xxx[2] contains the slip for second foc mech (for single patch events only!).
	double *slip_str;	//slip along strike
	double *slip_dip;	//slip along dip
//	double *strikes;	//strike of patch (not used)
//	double *dips;		//dip of patch (not used)
    double *pos_s;		//along strike distance from point: lat,lon (0 for single patch events);
    double *pos_d;		//along dip distance from point: lat,lon (0 for single patch events);
    double *distance;
    int *selpoints;		//indices of cell points affected by this event.
    int index_cat;		//index of event in catalog.
    struct set_of_models *parent_set_of_models;	//if multiple models are present, this is a pointer to corresponding set_of_models structure.
};

// [Fahad] A complete collection of scalar model parameters used in
//		   read_modelparmeters(). In terms of communications, it is
//		   more efficient to pack all the variables in one struct
//		   and then transport these to the other nodes over the network,
//		   rather than sending each variable separately.
struct BCast_Model_Parameters {
	int N_min_events;
	int fixr;
	int fixAsig;
	int fixta;
	int nAsig0;
	int nta0;
	int Nsur;
	int Nslipmod;
	int use_bg_rate;
	int gridPMax;
	int LLinversion;
	int forecast;
	double r0;
	double Asig0;
	double ta0;
	double Asig_min;
	double Asig_max;
	double ta_min;
	double ta_max;
	double tstartLL;
	double extra_time;
	double tw;
	double fore_dt;
	double t_back;
	double Hurst;
	double Mc_source;
	double Mc;
	double Mag_main;
	double DCFS_cap;
	double dt;
	double dM;
	double xytoll;
	double ztoll;
	double border;
	double res;
	double gridresxy;
	double gridresz;
	double smoothing;
};


#endif //DEFINES_H
