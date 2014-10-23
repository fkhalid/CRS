#!/bin/bash

if [ $# -eq 0 ]
then
 Build="Release"
else
 Build=$1
fi

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testC/input.txt"
parafile="input_testcases/parameters_uncert.txt"
temppara="input_testcases/testC/temp_par.txt"

for i in $(seq 0 4)	#should be 0 4
do
ln1="OutputForecastFile=output_testcases/testC$i"
ln2="Logfile=output_testcases/testC$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputC
sed "47s+0+1+"  $parafile | sed "47s+X+$i+"  > $temppara

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



