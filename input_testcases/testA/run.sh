#!/bin/bash

basefile="input_testcases/testA/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testA/temp_par.txt"

#run with vary_sm:
ln1="OutputForecastFile=output_testcases/testA1"
ln2="Logfile=output_testcases/testA1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
sed "41s+0+1+"  $parafile > $temppara

Release/CRS_3.0 temp_inputA

#run with vary_fm:
ln1="OutputForecastFile=output_testcases/testA2"
ln2="Logfile=output_testcases/testA2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
sed "42s+0+1+"  $parafile > $temppara

Release/CRS_3.0 temp_inputA

#run with vary_fm:
ln1="OutputForecastFile=output_testcases/testA3"
ln2="Logfile=output_testcases/testA3.log"
ln11="InputListCatalogFocMecFile=input_testcases/listfocmecfiles.txt"
sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "11s+.*+$ln11+" > temp_inputA
sed "42s+0+1+"  $parafile > $temppara

Release/CRS_3.0 temp_inputA

#run with vary_grid:
ln1="OutputForecastFile=output_testcases/testA4"
ln2="Logfile=output_testcases/testA4.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
sed "43s+0+1+"  $parafile > $temppara

Release/CRS_3.0 temp_inputA
rm temp_inputA
rm $temppara
