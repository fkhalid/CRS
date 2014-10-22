#!/bin/bash

if [ $# -eq 0 ]
then
 Build="Release"
else
 Build=$1
fi

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testG/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testG/temp_par.txt"

for i in $(seq 0 4)	#should be 0 4
do
if [ $i -eq 3 ]
then
continue
fi
ln1="OutputForecastFile=output_testcases/testG$i"
ln2="Logfile=output_testcases/testG$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputG
sed "47s+0+1+"  $parafile | sed "47s+X+$i+" | sed "67s+5.95+2.0+" > $temppara

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


#Here all events are treated as mainshocks, should be equivalent to testG1
for i in $(seq 0 4)     #should be 0 4
do
if [ $i -eq 3 ]
then
continue
fi
ln1="OutputForecastFile=output_testcases/testGB$i"
ln2="Logfile=output_testcases/testGB$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputG
sed "67s+.*+2.00+" $parafile | sed "47s+X+$i+" | sed "67s+5.95+2.0+" > $temppara

if [ $i -lt 2 ]
then
sed '9s+focmecfile+focmecfile0+' temp_inputG > tmp
mv tmp temp_inputG
fi


$Build/CRS_3.0 temp_inputG
if [ $Build == "Coverage" ]
then
mkdir coverage/testGB$i
cp Coverage/code/src/*/*.gc* coverage/testGB$i
fi

done

rm temp_inputG
rm $temppara



