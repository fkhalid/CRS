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
basefile="input_testcases/testDS/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testDS/temp_par.txt"

#no afterslip:
ln1="OutputForecastFile=output_testcases/testDS0"
ln2="Logfile=output_testcases/testDS0.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputDS
cp  $parafile $temppara
$Build/run_crs temp_inputDS
if [ $Build == "Coverage" ]
then
mkdir coverage/testDS0
cp Coverage/code/src/*/*.gc* coverage/testDS0
fi


#afterslip model with a single snapshot:
ln1="OutputForecastFile=output_testcases/testDS1"
ln2="Logfile=output_testcases/testDS1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputDS
echo "InputListAseismicModels=input_testcases/testDS/aslipmodels1.dat" >> temp_inputDS
cp  $parafile $temppara
$Build/run_crs temp_inputDS
if [ $Build == "Coverage" ]
then
mkdir coverage/testDS1
cp Coverage/code/src/*/*.gc* coverage/testDS1
fi


#afterslip model with multiple snapshots (splines):
ln1="OutputForecastFile=output_testcases/testDS2"
ln2="Logfile=output_testcases/testDS2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputDS
echo "InputListAseismicModels=input_testcases/testDS/aslipmodels2.dat" >> temp_inputDS
cp $parafile $temppara

$Build/run_crs temp_inputDS
if [ $Build == "Coverage" ]
then
mkdir coverage/testDS2
cp Coverage/code/src/*/*.gc* coverage/testDS2
fi


#afterslip model with a single snapshot:
ln1="OutputForecastFile=output_testcases/testDS3"
ln2="Logfile=output_testcases/testDS3.log"
ln13="InversionStartDate=2004-09-28T17:00:24Z"	#so it's comparable to testDS5

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "13s+.*+$ln13+" | sed "s*slipmodelslist0*slipmodelslist1*" > temp_inputDS
echo "InputListAseismicModels=input_testcases/testDS/aslipmodels3.dat" >> temp_inputDS
cp  $parafile $temppara
$Build/run_crs temp_inputDS
if [ $Build == "Coverage" ]
then
mkdir coverage/testDS3
cp Coverage/code/src/*/*.gc* coverage/testDS3
fi

#2 afterslips, with same geometry (single subfault)
ln1="OutputForecastFile=output_testcases/testDS4"
ln2="Logfile=output_testcases/testDS4.log"
ln7="InputCatalogFile=input_testcases/testDS/catalog2.dat"
ln13="InversionStartDate=2004-09-28T17:00:24Z"
ln6="ForecastEndDate=2005-10-15T17:15:24Z"
ln4="IssueDate=2005-09-30T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist2*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+" | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputDS
echo "InputListAseismicModels=input_testcases/testDS/aslipmodels4.dat" >> temp_inputDS
cp  $parafile $temppara
$Build/run_crs temp_inputDS
if [ $Build == "Coverage" ]
then
mkdir coverage/testDS4
cp Coverage/code/src/*/*.gc* coverage/testDS4
fi

#2 afterslip models, first one is empty (should be same as testDS3)
ln1="OutputForecastFile=output_testcases/testDS5"
ln2="Logfile=output_testcases/testDS5.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist1*" | sed "13s+.*+$ln13+" > temp_inputDS
echo "InputListAseismicModels=input_testcases/testDS/aslipmodels5.dat" >> temp_inputDS
cp  $parafile $temppara
$Build/run_crs temp_inputDS
if [ $Build == "Coverage" ]
then
mkdir coverage/testDS5
cp Coverage/code/src/*/*.gc* coverage/testDS5
fi

#2 afterslip models, with splines
ln1="OutputForecastFile=output_testcases/testDS6"
ln2="Logfile=output_testcases/testDS6.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist6*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+"  | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputDS
echo "InputListAseismicModels=input_testcases/testDS/aslipmodels6.dat" >> temp_inputDS
cp  $parafile $temppara
$Build/run_crs temp_inputDS
if [ $Build == "Coverage" ]
then
mkdir coverage/testDS6
cp Coverage/code/src/*/*.gc* coverage/testDS6
fi

