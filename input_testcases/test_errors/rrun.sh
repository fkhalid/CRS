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

echo "***************************build="$Build"********************************"
basefile="input_testcases/test_errors/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/test_errors/temp_par.txt"

#run with non chronological catalog:
ln1="OutputForecastFile=output_testcases/test_err1"
ln2="Logfile=output_testcases/test_err1.log"
ln7="InputCatalogFile=input_testcases/test_errors/catalog.dat"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "7s+.*+$ln7+"> temp_inputA
cp  $parafile  $temppara

$Build/run_crs temp_inputA
if [ $Build == "Coverage" ] 
then
mkdir coverage/test_err1
cp Coverage/code/src/*/*.gc* coverage/test_err1
fi

#run with non chronological slip model list:
ln1="OutputForecastFile=output_testcases/test_err2"
ln2="Logfile=output_testcases/test_err2.log"
ln8="InputListSlipModels=input_testcases/test_errors/slipmodelslist.dat"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "8s+.*+$ln8+"> temp_inputA
cp $parafile  $temppara

$Build/run_crs temp_inputA
if [ $Build == "Coverage" ]
then
mkdir coverage/test_err2
cp Coverage/code/src/*/*.gc* coverage/test_err2      
fi

#run with a line which is too long:
ln1="OutputForecastFile=output_testcases/test_err3_thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000thisisjustalinewhichisreallyreallylongandcannotbereadentirelyandshouldgiveanerrormessageaslongasNcharis1000"
ln2="Logfile=output_testcases/test_err3.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputA
cp $parafile  $temppara

$Build/run_crs temp_inputA
if [ $Build == "Coverage" ]
then
mkdir coverage/test_err3
cp Coverage/code/src/*/*.gc* coverage/test_err3
fi

#run with fixedmecfile file and wrong value for receiver fault orientation (it should be "fixed")
basefile="input_testcases/testF/input.txt"
parafile="input_testcases/testF/parameters_fixmec.txt"
temppara="input_testcases/testF/temp_par.txt"

ln1="OutputForecastFile=output_testcases/test_err3"
ln2="Logfile=output_testcases/test_err3.log"
sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputF
echo "FixedMecFile=input_testcases/testF/fixmec2.dat" >>  temp_inputF
echo "InputCatalogFocMecFile=input_testcases/focmecfile.dat" >> temp_inputF

sed "57s+.*+focmec+"  $parafile > $temppara

$Build/run_crs temp_inputF
if [ $Build == "Coverage" ]
then
mkdir coverage/testF3
cp Coverage/code/src/*/*.gc* coverage/testF3
fi


