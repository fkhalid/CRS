/*
 * mem_mgmt.h
 *
 *  Created on: Jun 5, 2013
 *      Author: camcat
 */

#ifndef MEM_MGMT_H_
#define MEM_MGMT_H_


#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "../defines.h"
#include "../util/nrutil.h"

void shift_cat(struct catalog *cat, int N);
void init_crst(struct crust *crst);
void init_cat1(struct catalog *cat, int Zsel);
struct set_of_models *set_of_models_array(long n1, long n2);
struct eqkfm *eqkfm_array(long n1, long n2);
struct pscmp *pscmp_array(long n1, long n2);
struct pscmp *pscmp_arrayinit(struct crust v0, long n1, long n2);
void free_eqkfmarray(struct eqkfm *v, long n1, long n2);
void freefull_eqkfmarray(struct eqkfm *v, long n1, long n2);
void freefull_pscmparray(struct pscmp *v, long n1, long n2);
void free_cat(struct catalog cat);

#endif /* MEM_MGMT_H_ */
