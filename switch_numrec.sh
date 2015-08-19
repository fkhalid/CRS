#!/bin/bash

old_numrec_folder="/home/des/camcat/masterdes_sec21/NumericalRecipes/2.11/blabla/"
new_numrec_folder="/home/des/camcat/masterdes_sec21/NumericalRecipes/2.11/C_211/"

grep $old_numrec_folder code/src/*/*mk -l | xargs sed -i "s*$old_numrec_folder*$new_numrec_folder*g"
