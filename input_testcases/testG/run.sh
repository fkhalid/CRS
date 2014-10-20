#!/bin/bash

basefile="input_testcases/testC/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testC/temp_par.txt"

for i in $(seq 0 4)	#should be 0 4
do
if [ $i -eq 3 ]
then
continue
fi
ln1="OutputForecastFile=output_testcases/testG$i"
ln2="Logfile=output_testcases/testG$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputC
sed "47s+0+1+"  $parafile | sed "47s+X+$i+" | sed "67s+5.95+2.0+" > $temppara

if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputC > tmp
mv tmp temp_inputC
fi

Release/CRS_3.0 temp_inputC
done


#Here all events are treated as mainshocks, should be equivalent to testC1
for i in $(seq 0 4)     #should be 0 4
do
if [ $i -eq 3 ]
then
continue
fi
ln1="OutputForecastFile=output_testcases/testGB$i"
ln2="Logfile=output_testcases/testGB$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputC
sed "67s+.*+2.00+" $parafile | sed "47s+X+$i+" | sed "67s+5.95+2.0+" > $temppara

if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputC > tmp
mv tmp temp_inputC
fi


Release/CRS_3.0 temp_inputC
done

rm temp_inputC
rm $temppara



