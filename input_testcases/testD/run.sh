#!/bin/bash

basefile="input_testcases/testD/input.txt"
parafile="input_testcases/parameters_aslip.txt"
temppara="input_testcases/testD/temp_par.txt"

#run with vary_sm:
ln1="OutputForecastFile=output_testcases/testD1"
ln2="Logfile=output_testcases/testD1.log"

#afterslip model with a single snapshot:
sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels1.dat" >> temp_inputD
sed "41s+0+1+"  $parafile > $temppara
Release/CRS_3.0 temp_inputD

#afterslip model with multiple snapshots (splines):
ln1="OutputForecastFile=output_testcases/testD2"
ln2="Logfile=output_testcases/testD2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels2.dat" >> temp_inputD
sed "42s+0+1+"  $parafile > $temppara

Release/CRS_3.0 temp_inputD
