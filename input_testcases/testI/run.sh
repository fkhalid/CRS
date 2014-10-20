#!/bin/bash

basefile="input_testcases/testI/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testI/temp_par.txt"

for i in $(seq 0 4)	#should be 0 4
do
if [ $i -eq 3 ]
then
continue
fi
ln1="OutputForecastFile=output_testcases/testI$i"
ln2="Logfile=output_testcases/testI$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputI
sed "47s+0+1+"  $parafile | sed "47s+X+$i+" | sed "67s+5.95+2.0+" > $temppara

if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputI > tmp
mv tmp temp_inputI
fi

Release/CRS_3.0 temp_inputI
done

rm temp_inputI
rm $temppara



