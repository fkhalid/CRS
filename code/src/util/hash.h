/*
 * hash.h
 *
 *  Created on: Sep 2, 2013
 *      Author: camcat
 */

#ifndef HASH_H_
#define HASH_H_


#endif /* HASH_H_ */


#include <stdio.h>      /* defines printf for tests */
#include <time.h>       /* defines time_t for timings in the test */
#include <stdint.h>     /* defines uint32_t etc */
#include <sys/param.h>  /* attempt to define endianness */
#ifdef linux
# include <endian.h>    /* attempt to define endianness */
#endif

uint32_t hashlittle( const void *key, size_t length, uint32_t initval);
