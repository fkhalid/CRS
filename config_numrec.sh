#/bin/bash

# This script sets the folder containing source code for Numerical Recipes.
#
# The old_numrec_folder is the current folder path. By default this is NumRec/folder/path. If it has been changed, you can find it by listing all include folders with this command:
# 	grep -o -r -h "\-I[1-9\.\/a-zA-Z_-]*" Release/ | sort | uniq
# for example, if the output of this command is:
# -I/NR/folder/other
# -I/NR/folder/recipes
# you should set old_numrec_folder="NR/folder/"

# The new_numrec_folder is the folder containing the Numerical Recipes source code. This folder should contain a "/recipes" and an "/others" subfolder.
# The format of this variable should be "/path/to/numericalrecipes/"

old_numrec_folder="/path/to/NR/"
new_numrec_folder="/new/path/to/NR/"

if [ $(grep -rl "$old_numrec_folder" Release | wc -l) == 0 ];
 then
        echo
	echo " The original path folder for numerical recipes ($old_numrec_folder) was not found in the makefiles."
	echo " Please set correct path in config_numrec.sh."
	echo " The command 'grep -o -r -h \"\-I[1-9\.\/a-zA-Z_-]*\" Release/ | sort | uniq' can be used to list the current include paths."
	echo
	echo " Numerical recipes path was not changed.";
	echo
        exit
 else

 # This commands finds all occurrences of the old include and overwrites it with the correct one:
 grep $old_numrec_folder Release -lr | xargs sed -i "s*$old_numrec_folder*$new_numrec_folder*g"

 echo 
 echo " The following include patterns have been set:"
 grep -o -r -h "\-I[1-9\.\/a-zA-Z_-]*" Release/ | sort | uniq
 echo

fi





