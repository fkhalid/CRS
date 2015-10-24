#!/bin/bash

grep "\.\.\/util\/ran1\.h" code/src/*/* -l | xargs sed -i "s*\.\.\/util\/ran1\.h*nr\.h*g"
grep "\.\.\/util\/nr\.h" code/src/*/* -l | xargs sed -i "s*\.\.\/util\/nr\.h*nr\.h*g"
grep "\.\.\/util\/nrutil\.h" code/src/*/* -l | xargs sed -i "s*\.\.\/util\/nrutil\.h*nrutil\.h*g"
