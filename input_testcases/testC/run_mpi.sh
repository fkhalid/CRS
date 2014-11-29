#!/bin/bash

Build="mpirun -n 2 Release/"

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testC/input.txt"
parafile="input_testcases/parameters_uncert.txt"
temppara="input_testcases/testC/temp_par.txt"

m0="iso iso"
m1="fm  iso"
m2="fm fix"
m3="fm no"	
for i in $(seq 0 3)	#should be 0 3
do
ln1="OutputForecastFile=output_testcases/testC$i"
ln2="Logfile=output_testcases/testC$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputC
l160=`echo m$i`
l16=`echo ${!l160}`
sed "16s+.*+$l16+"  $parafile | sed '15s+5.95+2.00+'> $temppara
# use foc mec file which contains one of the source events (to distinguish test0 and test1).
#if [ $i -lt 2 ]
if [ $i -eq 1 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputC > tmp
mv tmp temp_inputC
fi

$Build/CRS_3.0 temp_inputC
if [ $Build == "Coverage" ]
then
mkdir coverage/testC$i
cp Coverage/code/src/*/*.gc* coverage/testC$i
fi

done

rm temp_inputC
rm $temppara



