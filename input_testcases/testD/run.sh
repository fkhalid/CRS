#!/bin/bash

if [ $# -eq 0 ]
then
 Build="Release"
else
 Build=$1
fi

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testD/input.txt"
parafile="input_testcases/parameters_aslip.txt"
temppara="input_testcases/testD/temp_par.txt"

#no afterslip:
ln1="OutputForecastFile=output_testcases/testD0"
ln2="Logfile=output_testcases/testD0.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputD
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD0
cp Coverage/code/src/*/*.gc* coverage/testD0
fi


#afterslip model with a single snapshot:
ln1="OutputForecastFile=output_testcases/testD1"
ln2="Logfile=output_testcases/testD1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels1.dat" >> temp_inputD
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD1
cp Coverage/code/src/*/*.gc* coverage/testD1
fi


#afterslip model with multiple snapshots (splines):
ln1="OutputForecastFile=output_testcases/testD2"
ln2="Logfile=output_testcases/testD2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels2.dat" >> temp_inputD
cp $parafile $temppara

$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD2
cp Coverage/code/src/*/*.gc* coverage/testD2
fi


#afterslip model with a single snapshot:
ln1="OutputForecastFile=output_testcases/testD3"
ln2="Logfile=output_testcases/testD3.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist1*" > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels3.dat" >> temp_inputD
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD3
cp Coverage/code/src/*/*.gc* coverage/testD3
fi

