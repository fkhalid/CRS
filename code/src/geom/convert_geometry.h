
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


#ifndef CONVERT_GEOMETRY_H_
#define CONVERT_GEOMETRY_H_

#include <stdio.h>

#include "../defines.h"

#include "../util/nrutil_newnames.h"

int convert_geometry(struct crust crst, double *old_v, double **new_v, int sum, int increase_resolution);
int flatten_outgrid(struct crust crst, double *old_v, double **new_v, int *Nfinal);

#endif /* CONVERT_GEOMETRY_H_ */
