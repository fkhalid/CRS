#!/bin/bash

basefile="input_testcases/testH/input.txt"
parafile="input_testcases/testH/parameters_as.txt"	#can use parameters_(m)as well to treat events as (main)aftershocks.
temppara="input_testcases/testH/temp_par.txt"

# considering events as mainshocks:
for i in $(seq 0 4)	#should be 0 4
do
ln1="OutputForecastFile=output_testcases/testH$i"
ln2="Logfile=output_testcases/testH$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputH
sed "47s+X+$i+" $parafile > $temppara

# use foc mec file which contains one of the source events (to distinguish test0 and test1).
if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputH > tmp
mv tmp temp_inputH
fi

Release/CRS_3.0 temp_inputH
done
rm temp_inputH
rm $temppara



