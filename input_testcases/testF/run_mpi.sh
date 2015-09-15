#!/bin/bash

Build="mpirun -n 4 Release/run_mpi_crs"

echo "***************************Build="$Build"********************************"
basefile="input_testcases/testF/input.txt"
parafile="input_testcases/testF/parameters_fixmec.txt"
temppara="input_testcases/testF/temp_par.txt"

#run without fixedmecfile file:
ln1="OutputForecastFile=output_testcases/testF1"
ln2="Logfile=output_testcases/testF1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputF
cp  $parafile $temppara
echo "InputCatalogFocMecFile=input_testcases/focmecfile.dat" >>  temp_inputF

$Build temp_inputF
if [ $Build == "Coverage" ]
then
mkdir coverage/testF1
cp Coverage/code/src/*/*.gc* coverage/testF1
fi


#run with fixedmecfile file with same mec. as above:
ln1="OutputForecastFile=output_testcases/testF2"
ln2="Logfile=output_testcases/testF2.log"
sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputF
cp  $parafile $temppara
echo "FixedMecFile=input_testcases/testF/fixmec1.dat" >>  temp_inputF

$Build temp_inputF
if [ $Build == "Coverage" ]
then
mkdir coverage/testF2
cp Coverage/code/src/*/*.gc* coverage/testF2
fi

#run with fixedmecfile file with random foc. planes:
ln1="OutputForecastFile=output_testcases/testF3"
ln2="Logfile=output_testcases/testF3.log"
sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputF
cp  $parafile $temppara
echo "FixedMecFile=input_testcases/testF/fixmec2.dat" >>  temp_inputF

$Build temp_inputF
if [ $Build == "Coverage" ]
then
mkdir coverage/testF3
cp Coverage/code/src/*/*.gc* coverage/testF3
fi

#run with fixedmecfile file with random foc. planes and refined grid:
#slow! only run once in a while
exit
ln1="OutputForecastFile=output_testcases/testF4"
ln2="Logfile=output_testcases/testF4.log"
sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputF
sed "17s+.*+2.0+" $parafile | sed "19s+.*+2.0    3.0+" > $temppara
echo "FixedMecFile=input_testcases/testF/fixmec2.dat" >>  temp_inputF

$Build temp_inputF
if [ $Build == "Coverage" ]
then
mkdir coverage/testF4
cp Coverage/code/src/*/*.gc* coverage/testF4
fi

rm $temppara

