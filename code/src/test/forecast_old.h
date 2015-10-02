
/*   Copyright (C) 2015 by Camilla Cattania and Fahad Khalid.
 *
 *   This file is part of CRS.
 *
 *   CRS is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CRS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CRS.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 * forecast_stepG.h
 *
 *  Created on: Jul 20, 2012
 *      Author: camcat
 */
#ifndef FORECAST_STEPG_H_
#define FORECAST_STEPG_H_

#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "../defines.h"

#include "../util/nrutil_newnames.h"
#include "../general/struct_conversions.h"

//void forecast_stepG(double *times, double **cmpdata, double tt0, double tt1, double Asig, double ta, int points[], double *NeX, double *NeT, int N, int NTS, double *gamma_init, int last);
//void forecast_stepG2(double *times, double **cmpdata, double tt0, double tt1, double Asig, double ta, int points[], double *NeX, double *NeT, int N, int NTS, double *gamma_init, int last);
int forecast_stepG2_old(struct catalog cat, double *times, double **cmpdata, struct pscmp *DCFS, double tt0, double tt1, double Asig, double ta, int points[], double *NeX, double *NeT, double *Rate_end, int N, int NTS, int Neq, double *gamma_init, double *back_rate, double *R, int last);
//void forecast_stepG2_iso(struct catalog cat, double *times, double **cmpdata, struct dist *DIST, struct ETASparameter ETAS, double tt0, double tt1, double Asig, double ta, int points[], double *NeX, double *NeT, int N, int NTS, int Neq, double *gamma_init, double *R, int last);

#endif /* FORECAST_STEPG_H_ */
