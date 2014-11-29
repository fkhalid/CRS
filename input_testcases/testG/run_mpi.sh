#!/bin/bash

Build="mpirun -n 2 Release/"

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testG/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testG/temp_par.txt"

m0="iso iso"
m1="fm  iso"
m2="fm fix"
m3="fm no"
for i in $(seq 0 3)     #should be 0 3
do
ln1="OutputForecastFile=output_testcases/testG$i"
ln2="Logfile=output_testcases/testG$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputG
l160=`echo m$i`
l16=`echo ${!l160}`
sed "16s+.*+$l16+"  $parafile | sed '15s+5.95+2.00+'> $temppara

if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputG > tmp
mv tmp temp_inputG
fi

$Build/CRS_3.0 temp_inputG
if [ $Build == "Coverage" ]
then
mkdir coverage/testG$i
cp Coverage/code/src/*/*.gc* coverage/testG$i
fi

done


#Here events are tapered, should be similar to G1

ln1="OutputForecastFile=output_testcases/testGT1"
ln2="Logfile=output_testcases/testGT1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputG
sed "16s+.*+$m1+"  $parafile | sed '15s+5.95+2.00+' | sed "17s+.*+3.0+" > $temppara
sed '9s+focmecfile+focmecfile0+' temp_inputG > tmp
mv tmp temp_inputG

$Build/CRS_3.0 temp_inputG
if [ $Build == "Coverage" ]
then
mkdir coverage/testGT1
cp Coverage/code/src/*/*.gc* coverage/testGT1
fi

rm temp_inputG
rm $temppara



