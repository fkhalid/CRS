/*
 * eqkfm_copy.h
 *
 *  Created on: Nov 20, 2013
 *      Author: camcat
 */

#ifndef EQKFM_COPY_H_
#define EQKFM_COPY_H_

#include "../defines.h"

void empty_eqkfm(struct eqkfm *eqkfm0);
void copy_eqkfm_nolocation_noindex_notime(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_attributes(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_focmec(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_all(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);
void copy_eqkfm_noslipmodel(struct eqkfm eqkfm1, struct eqkfm *eqkfm2);

#endif /* EQKFM_COPY_H_ */
