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
basefile="input_testcases/testD/input.txt"
parafile="input_testcases/parameters.txt"
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
ln13="InversionStartDate=2004-09-28T17:00:24Z"	#so it's comparable to testD5

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "13s+.*+$ln13+" | sed "s*slipmodelslist0*slipmodelslist1*" > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels3.dat" >> temp_inputD
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD3
cp Coverage/code/src/*/*.gc* coverage/testD3
fi

#2 afterslips, with same geometry (single subfault)
ln1="OutputForecastFile=output_testcases/testD4"
ln2="Logfile=output_testcases/testD4.log"
ln7="InputCatalogFile=input_testcases/testD/catalog2.dat"
ln13="InversionStartDate=2004-09-28T17:00:24Z"
ln6="ForecastEndDate=2005-10-15T17:15:24Z"
ln4="IssueDate=2005-09-30T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist2*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+" | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels4.dat" >> temp_inputD
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD4
cp Coverage/code/src/*/*.gc* coverage/testD4
fi

#2 afterslip models, first one is empty (should be same as testD3)
ln1="OutputForecastFile=output_testcases/testD5"
ln2="Logfile=output_testcases/testD5.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist1*" | sed "13s+.*+$ln13+" > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels5.dat" >> temp_inputD
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD5
cp Coverage/code/src/*/*.gc* coverage/testD5
fi

#2 afterslip models, with splines
ln1="OutputForecastFile=output_testcases/testD6"
ln2="Logfile=output_testcases/testD6.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist2*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+"  | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputD
echo "InputListAfterslipModels=input_testcases/testD/aslipmodels6.dat" >> temp_inputD
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputD
if [ $Build == "Coverage" ]
then
mkdir coverage/testD6
cp Coverage/code/src/*/*.gc* coverage/testD6
fi

