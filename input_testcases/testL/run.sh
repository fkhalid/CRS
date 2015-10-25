#!/bin/bash

if [ $# -eq 0 ]
then
 Build="Release"
else
 if [ $1 == "MPI" ]
 then
  Build="mpirun -n $2 MPI/"
 else
  Build=$1
 fi
fi

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testL/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testL/temp_par.txt"

for i in $(seq 1 4)
do
ln1="OutputForecastFile=output_testcases/testL$i"
ln2="Logfile=output_testcases/testL$i.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputL
echo "InputListSlipModels=input_testcases/testL/slipmodelslist$i.dat" >> temp_inputL
cp $parafile $temppara

$Build/run_crs temp_inputL
if [ $Build == "Coverage" ]
then
mkdir coverage/testL$i
cp Coverage/code/src/*/*.gc* coverage/testL$i
fi
done

rm temp_inputL
rm $temppara

