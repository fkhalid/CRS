#/bin/bash

# This script sets the folder containing source code for Numerical Recipes, or allows to compile the code without using them.
# You only need the Numerical Recipes if you want to use the option 'splines' for aseismic stress sources. If you are not including aseismic processes, or if you use another option ('log' or 'lin' in the parameter file), you can set: "numerical_recipes_installed=no" and run this script. Otherwise, you need to set the paths to the Numerical Recipes source files below.
# We hope to make the code completely Numerical-Recipes free in the future!
#
# The old_numrec_folder is the current folder path, which is found at the end of Release/code/src/*/subdir.mk files. By default this is NumRec/folder/path/. If it has been changed, you can find it by listing all include folders with this command:
# 	grep -o -r -h "\-I[1-9\.\/a-zA-Z_-]*" Release/ | sort | uniq
# for example, if the output of this command is:
# -I/NR/folder/other
# -I/NR/folder/recipes
# you should set old_numrec_folder="NR/folder/"
# The format of this variable should be "/path/to/numericalrecipes/"


# Write yes or no here:
numerical_recipes_installed=no

# Put the old path to the Numerical Recipes source code here:
old_numrec_folder="/home/des/camcat/masterdes_exthd1/sec21/NumericalRecipes/2.11/C_211/"
# Put the path to the Numerical Recipes source code here:
new_numrec_folder="/home/des/camcat/masterdes_exthd1/sec21/NumericalRecipes/2.11/C_212/"


#---------------Don't need to change anything below here--------------------------------%


# Name of flag to deactivae numerical recipes.
noNRflag="_no_numericalrecipes"

if [ $numerical_recipes_installed = "yes" ];
  then
	#Remove flag which deactivates numerical recipes:
        grep -rl "$noNRflag" Release | xargs sed -i "s* -D$noNRflag**g"
        echo 
	echo " Removed flag $noNRflag from makefiles."
	echo

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
elif [ $numerical_recipes_installed = "no" ];
  then

	#determine which compiler is set (gcc or mpicc)

	find_gcc=$(grep -rl gcc Release/ | grep 'subdir\|makefile' | wc -l)
	find_mpicc=$(grep -rl mpicc Release/ | grep 'subdir\|makefile' | wc -l)

	if [ $find_gcc != "0" ] && [ $find_mpicc != "0" ]
	then
	   echo 
	   echo " Could not determine which compiler is being used: please set it manually in config_numrec.sh"
	   echo
	else
	   if [ $find_gcc != "0" ]
	   then
		grep -rl gcc Release/ | grep 'subdir\.mk' | xargs grep -rL $noNRflag | xargs sed -i "s*gcc*gcc -D$noNRflag*"
	   else
		grep -rl mpicc Release/ | grep 'subdir\.mk' | xargs grep -rL $noNRflag | xargs sed -i "s*mpicc*mpicc -D$noNRflag*"
	   fi
	   echo
	   echo " Added flag: $noNRflag to makefiles."
	   echo " Now you can compile without the Numerical Recipes source code."
	   echo 	
	fi	

else
  echo
  echo " Illegal value for variable 'numerical_recipes_installed' in config_numrec.sh. Should be 'yes' or 'no' "
  echo
fi   
  
  
  
  
