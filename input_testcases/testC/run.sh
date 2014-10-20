#!/bin/bash

basefile="input_testcases/testC/input.txt"
parafile="input_testcases/parameters_uncert.txt"
temppara="input_testcases/testC/temp_par.txt"

for i in $(seq 0 4)	#should be 0 4
do
ln1="OutputForecastFile=output_testcases/testC$i"
ln2="Logfile=output_testcases/testC$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputC
sed "47s+0+1+"  $parafile | sed "47s+X+$i+"  | sed "67s+5.95+2.0+" > $temppara

# use foc mec file which contains one of the source events (to distinguish test0 and test1).
#if [ $i -lt 2 ]
if [ $i -eq 1 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputC > tmp
mv tmp temp_inputC
fi

Release/CRS_3.0 temp_inputC
done

rm temp_inputC
rm $temppara



