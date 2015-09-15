#!/bin/bash

Build="mpirun -n 4 Release/run_mpi_crs"

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testB/input.txt"
parafile="input_testcases/parameters_uncert.txt"
temppara="input_testcases/testB/temp_par.txt"

#run with catalog2 (mainshock at t=0 moved back by 100 days):
ln1="OutputForecastFile=output_testcases/testB1"
ln2="Logfile=output_testcases/testB1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "7s+catalog+catalog2+" > temp_inputB
cp  $parafile $temppara

$Build temp_inputB
if [ $Build == "Coverage" ]
then
mkdir coverage/testB1
cp Coverage/code/src/*/*.gc* coverage/testB1
fi

#run with catalog2 (mainshock at t=0 removed):
ln1="OutputForecastFile=output_testcases/testB2"
ln2="Logfile=output_testcases/testB2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "7s+catalog+catalog3+" > temp_inputB
cp $parafile $temppara

$Build temp_inputB
if [ $Build == "Coverage" ]
then
mkdir coverage/testB2
cp Coverage/code/src/*/*.gc* coverage/testB2
fi

rm temp_inputB
rm $temppara
