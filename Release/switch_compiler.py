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
    selfFile = './' + makeFile
    print 'Updating: ' + makeFile
    
    f = open(makeFile,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(sys.argv[1],sys.argv[2])

    f = open(makeFile,'w')
    f.write(newdata)
    f.close()
