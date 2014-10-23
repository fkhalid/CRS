#!/bin/bash

if [ $# -eq 0 ]
then
 Build="Release"
else
 Build=$1
fi

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testM/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testM/temp_par.txt"

for i in $(seq 1 2)
do
ln1="OutputForecastFile=output_testcases/testM$i"
ln2="Logfile=output_testcases/testM$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputM
echo "ForecastTemplate=input_testcases/darf_temp_fake$i.txt" >> temp_inputM
cp $parafile $temppara

$Build/CRS_3.0 temp_inputM
if [ $Build == "Coverage" ]
then
mkdir coverage/testM$i
cp Coverage/code/src/*/*.gc* coverage/testM$i
fi
done

rm temp_inputM
rm $temppara

