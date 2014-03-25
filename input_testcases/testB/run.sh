#!/bin/bash

basefile="input_testcases/testB/input.txt"
parafile="input_testcases/parameters_uncert.txt"
temppara="input_testcases/testB/temp_par.txt"

#run with catalog2 (mainshock at t=0 moved back by 100 days):
ln1="OutputForecastFile=output_testcases/testB1"
ln2="Logfile=output_testcases/testB1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "7s+catalog+catalog2+" > temp_inputB
cp  $parafile $temppara

Release/CRS_3.0 temp_inputB

#run with catalog2 (mainshock at t=0 removed):
ln1="OutputForecastFile=output_testcases/testB2"
ln2="Logfile=output_testcases/testB2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "7s+catalog+catalog3+" > temp_inputB
cp $parafile $temppara

Release/CRS_3.0 temp_inputB

rm temp_inputB
rm $temppara
