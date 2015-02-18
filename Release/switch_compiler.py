#!/usr/bin/env python

import os
import sys
import fileinput

makeFiles = os.popen("grep -rl " + sys.argv[1] + " *").read().split()
#makeFile = os.popen("grep -rl " + sys.argv[1]).read().split()	//"grep -lr <blabla>" (without the *) hangs...


if len(makeFiles) == 0 :
    print 'No file found with string: ' + sys.argv[1]
else :
    print str(len(makeFiles)) + ' files found with string: ' + sys.argv[1]


for makeFile in makeFiles :    
    f = open(makeFile,'r')
    filedata = f.read()
    f.close()

    tokens = makeFile.split('/');

    if sys.argv[2] == "mpicc" :
        # Set the _CRS_MPI directive when compiling for MPI.
        if "subdir.mk" in tokens :
            # All subdir.mk files are used for compiling, and not linking.
            # Therefore, the -D_CRS_MPI flag should be added.
            newdata = filedata.replace(sys.argv[1], sys.argv[2] + " -D_CRS_MPI")
        else :
            # makefile is used for linking only. Therefore, the -D_CRS_MPI
            # should not be added.
            newdata = filedata.replace(sys.argv[1], sys.argv[2])
    else :
        # The _CRS_MPI should not be set when compiling serial code.
	    if "subdir.mk" in tokens :
            	newdata = filedata.replace(sys.argv[1] + " -D_CRS_MPI", sys.argv[2])
            else :
                newdata = filedata.replace(sys.argv[1], sys.argv[2])

    # Just to make sure this script itself is not updated.
    if "switch_compiler.py" not in tokens :
        f = open(makeFile,'w')
        f.write(newdata)
        f.close()
        print 'Updated: ' + makeFile
