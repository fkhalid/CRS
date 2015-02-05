#!/bin/bash

Build="mpirun -n 2 MPI/"

echo "***************************build="$Build"********************************"
basefile="input_testcases/testA/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testA/temp_par.txt"

#run with vary_sm:
#ln1="OutputForecastFile=output_testcases/testA1"
#ln2="Logfile=output_testcases/testA1.log"

#sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
#sed "43s+0+1+"  $parafile > $temppara
#echo "InputCatalogFocMecFile=input_testcases/focmecfile.dat" >> temp_inputA
#$Build/CRS_3.0 temp_inputA

#run with vary_fm:
ln1="OutputForecastFile=output_testcases/testA2"
ln2="Logfile=output_testcases/testA2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
#sed "43s+0+1+"  $parafile > $temppara
sed "57s+.*+focmec+" $parafile > $temppara
echo "InputCatalogFocMecFile=input_testcases/focmecfile.dat" >> temp_inputA

$Build/CRS_3.0 temp_inputA
if [ $Build == "Coverage" ] 
then
mkdir coverage/testA2
cp Coverage/code/src/*/*.gc* coverage/testA2
fi

#run with vary_fm:
ln1="OutputForecastFile=output_testcases/testA3"
ln2="Logfile=output_testcases/testA3.log"
sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
sed "57s+.*+focmec+" $parafile > $temppara
#sed "43s+0+1+"  $parafile > $temppara
echo "InputListCatalogFocMecFile=input_testcases/listfocmecfiles.txt" >> temp_inputA

$Build/CRS_3.0 temp_inputA
if [ $Build == "Coverage" ]
then
mkdir coverage/testA3
cp Coverage/code/src/*/*.gc* coverage/testA3      
fi

#run with vary_grid:
ln1="OutputForecastFile=output_testcases/testA4"
ln2="Logfile=output_testcases/testA4.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
#sed "44s+0+1+"  $parafile > $temppara
sed "58s+.*+1+" $parafile > $temppara
echo "InputCatalogFocMecFile=input_testcases/focmecfile.dat" >> temp_inputA

$Build/CRS_3.0 temp_inputA
if [ $Build == "Coverage" ]
then
mkdir coverage/testA4
cp Coverage/code/src/*/*.gc* coverage/testA4
fi


rm temp_inputA
rm $temppara
