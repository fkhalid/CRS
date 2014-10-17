#!/bin/bash

basefile="input_testcases/testC/input.txt"
parafile="input_testcases/parameters_uncert.txt"
temppara="input_testcases/testC/temp_par.txt"

for i in $(seq 0 4)	#should be 0 4
do
ln1="OutputForecastFile=output_testcases/testC$i"
ln2="Logfile=output_testcases/testC$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputC
sed "49s+0+1+"  $parafile | sed "49s+X+$i+" > $temppara
mpirun -np 2 Release/CRS_3.0 temp_inputC

rm temp_inputC
rm $temppara
done


