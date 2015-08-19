#!/bin/bash

folder=Release;

if [ $1 == "gcc2mpicc" ]
then

 old_compiler=gcc;
 new_compiler=mpicc;
 old_name=run_crs;
 new_name=run_crs_mpi;

 # Output message and do nothing if compiler not found:
 if [ $(grep -rl "$old_compiler" $folder | grep -v -x ".*\.py\|.*\.sh" | wc -l) == 0 ];
 then 
	echo "No files with compiler $old_compiler found.";
	exit
 else

   # All subdir.mk files are used for compiling, and not linking.
   # Therefore, the -D_CRS_MPI flag should be added.
   grep -rl "$old_compiler" $folder | grep "subdir.mk" | xargs sed -i "s*$old_compiler*$new_compiler -D_CRS_MPI*"
  
   # makefile is used for linking only. Therefore, the -D_CRS_MPI
   # should not be added.
   #  grep -v -x ".*\.py\|.*\.sh" to avoid changing variables in scripts.
   grep -rl "$old_compiler" $folder | grep -v "subdir.mk" |  grep -v -x ".*\.py\|.*\.sh" | xargs sed -i "s*$old_compiler*$new_compiler*"

   #change name of the executable:
   grep -rl "$old_name" $folder |  grep -v -x ".*\.py\|.*\.sh" | xargs sed -i "s*$old_name*$new_name*"

   #output message
   echo "Switched compiler gcc with mpicc."
 fi

else
 if [ $1 == "mpicc2gcc" ]
 then

 old_compiler=mpicc;
 new_compiler=gcc; 
 old_name=run_crs_mpi;
 new_name=run_crs;


 if [ $(grep -rl "$old_compiler" $folder | grep -v -x ".*\.py\|.*\.sh" | wc -l) == 0 ];
 then  
        echo "No files with compiler $old_compiler found.";
        exit
 else 

   # The _CRS_MPI should not be set when compiling serial code. 
   grep -rl "$old_compiler" $folder | grep "subdir.mk" | xargs sed -i "s*$old_compiler -D_CRS_MPI*$new_compiler*"
   grep -rl "$old_compiler" $folder | grep -v "subdir.mk" |  grep -v -x ".*\.py\|.*\.sh" | xargs sed -i "s*$old_compiler*$new_compiler*"

   echo "Switched compiler mpicc with gcc."
 fi
 
 else
   echo "Usage:";
   echo "   switch_compiler gcc2mpicc (to go from serial to parallel)"
   echo "   switch_compiler mpicc2gcc (to go from parallel to serial)"
 fi
fi

 

