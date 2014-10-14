#!/bin/bash

basefile="input_testcases/testE/input.txt"
parafile="input_testcases/parameters_uncert.txt"
temppara="input_testcases/testE/temp_par.txt"

#run with internal background rate calculations (but not enough events available):
ln1="OutputForecastFile=output_testcases/testE1"
ln2="Logfile=output_testcases/testE1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputE
sed "48s+0+1+"  $parafile > $temppara
echo "InputBackgroundRateCatalog=input_testcases/catalog.dat" >> temp_inputE
Release/CRS_3.0 temp_inputE

#run with internal background rate calculations (not enough events available after declustering):
ln1="OutputForecastFile=output_testcases/testE2"
ln2="Logfile=output_testcases/testE2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "7s+catalog+catalog4+" > temp_inputE
sed "48s+0+1+"  $parafile > $temppara
echo "InputBackgroundRateCatalog=input_testcases/catalog4.dat" >> temp_inputE
Release/CRS_3.0 temp_inputE

#run with internal background rate calculations (enough events available):
ln1="OutputForecastFile=output_testcases/testE3"
ln2="Logfile=output_testcases/testE3.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "7s+catalog+catalog5b+" > temp_inputE
sed "48s+0+1+"  $parafile > $temppara
echo "InputBackgroundRateCatalog=input_testcases/catalog5b.dat" >> temp_inputE

Release/CRS_3.0 temp_inputE


#run with externally calculated background rate (same as previous case):
ln1="OutputForecastFile=output_testcases/testE4"
ln2="Logfile=output_testcases/testE4.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "7s+catalog+catalog5b+" > temp_inputE
echo "InputBackgroundRateGrid=input_testcases/background3.dat" >> temp_inputE

Release/CRS_3.0 temp_inputE

rm temp_inputE
rm $temppara
