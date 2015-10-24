#!/bin/bash

#grep -rl nrutil_newnames code | xargs sed -i 's*#include "../util/nrutil_newnames.h"*#include "../util/nrutil_newnames.h"\n#include "nrutil.h"*g'
#grep -rl '#include "nrutil.h"' code | xargs sed -i 's*#include "nrutil.h"*#include "../util/nrutil.h"*g'

grep -rl dvector code | grep -v "nrutil" | xargs sed -i 's*dvector*darray*g'
#grep -rl vector code| grep -v "nrutil"  | xargs sed -i 's*vector*farray*g'
grep -rl ivector code | grep -v "nrutil" | xargs sed -i 's*ivector*iarray*g'
#grep -rl larray code | grep -v "newnames" | xargs sed -i 's*larray*lvector*g'
#grep -rl f3array code | grep -v "newnames" | xargs sed -i 's*f3array*f3tensor*g'
#grep -rl d2array code | grep -v "newnames" | xargs sed -i 's*d2array*dmatrix*g'
#grep -rl f2array code | grep -v "newnames" | xargs sed -i 's*f2array*matrix*g'
#grep -rl i2array code | grep -v "newnames" | xargs sed -i 's*i2array*imatrix*g'
