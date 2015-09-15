#!/bin/bash

Build="mpirun -n 4 Release/run_mpi_crs"

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testH/input.txt"
parafile="input_testcases/testH/parameters_as.txt"	#can use parameters_(m)as well to treat events as (main)aftershocks.
temppara="input_testcases/testH/temp_par.txt"

m0="iso iso"
m1="fm  iso"
m2="fm fix"
m3="fm no"
for i in $(seq 0 3)     #should be 0 3
do
ln1="OutputForecastFile=output_testcases/testH$i"
ln2="Logfile=output_testcases/testH$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputH
l160=`echo m$i`
l16=`echo ${!l160}`
sed "16s+.*+$l16+"  $parafile | sed '15s+5.95+2.00+'> $temppara

# use foc mec file which contains one of the source events (to distinguish test0 and test1).
if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputH > tmp
mv tmp temp_inputH
fi

$Build temp_inputH
if [ $Build == "Coverage" ]
then
mkdir coverage/testH$i
cp Coverage/code/src/*/*.gc* coverage/testH$i
fi

done
rm temp_inputH
rm $temppara



