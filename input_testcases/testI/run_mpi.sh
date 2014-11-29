#!/bin/bash

Build="mpirun -n 2 Release/"

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testI/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testI/temp_par.txt"

m0="iso iso"
m1="fm  iso"
m2="fm fix"
m3="fm no"
for i in $(seq 0 3)     #should be 0 3
do
ln1="OutputForecastFile=output_testcases/testI$i"
ln2="Logfile=output_testcases/testI$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputI
l160=`echo m$i`
l16=`echo ${!l160}`
sed "16s+.*+$l16+"  $parafile | sed '15s+5.95+2.00+'> $temppara

#use a different input file with foc. mec. for the mainshocks:
if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputI  > tmp
mv tmp temp_inputI
fi
$Build/CRS_3.0 temp_inputI
if [ $Build == "Coverage" ]
then
mkdir coverage/testI$i
cp Coverage/code/src/*/*.gc* coverage/testI$i
fi
done

#Test a model in which the three slip models are identical, and one with a single slip model (should give same result):

ln1="OutputForecastFile=output_testcases/testIB1"
ln2="Logfile=output_testcases/testIB1.log"
ln8="InputListSlipModels=input_testcases/slipmodelslist.dat"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "8s+.*+$ln8+" > temp_inputI
sed "16s+.*+$m0+"  $parafile | sed '15s+5.95+2.00+'> $temppara

$Build/CRS_3.0 temp_inputI
if [ $Build == "Coverage" ]
then
mkdir coverage/testIB1
cp Coverage/code/src/*/*.gc* coverage/testIB1
fi

#3 slip models with same geometry:
ln1="OutputForecastFile=output_testcases/testIB2"
ln2="Logfile=output_testcases/testIB2.log"
ln8="InputListSlipModels=input_testcases/testI/slipmodelslist2.dat"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "8s+.*+$ln8+" > temp_inputI
sed "16s+.*+$m0+"  $parafile | sed '15s+5.95+2.00+'> $temppara

$Build/CRS_3.0 temp_inputI
if [ $Build == "Coverage" ]
then
mkdir coverage/testIB2
cp Coverage/code/src/*/*.gc* coverage/testIB2
fi

#here pretend slip  models have different geometry:
ln1="OutputForecastFile=output_testcases/testIB3"
ln2="Logfile=output_testcases/testIB3.log"
ln8="InputListSlipModels=input_testcases/testI/slipmodelslist3.dat"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "8s+.*+$ln8+" > temp_inputI
sed "16s+.*+$m0+"  $parafile | sed '15s+5.95+2.00+'> $temppara

$Build/CRS_3.0 temp_inputI
if [ $Build == "Coverage" ]
then
mkdir coverage/testIB3
cp Coverage/code/src/*/*.gc* coverage/testIB3
fi


rm temp_inputI
rm $temppara

