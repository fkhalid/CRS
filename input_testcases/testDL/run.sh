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
basefile="input_testcases/testDL/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testDL/temp_par.txt"

#no afterslip:
ln1="OutputForecastFile=output_testcases/testDL0"
ln2="Logfile=output_testcases/testDL0.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL0
cp Coverage/code/src/*/*.gc* coverage/testDL0
fi


#afterslip model with a single snapshot:
ln1="OutputForecastFile=output_testcases/testDL1"
ln2="Logfile=output_testcases/testDL1.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels1.dat" >> temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL1
cp Coverage/code/src/*/*.gc* coverage/testDL1
fi


#afterslip model with multiple snapshots (splines):
ln1="OutputForecastFile=output_testcases/testDL2"
ln2="Logfile=output_testcases/testDL2.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels2.dat" >> temp_inputDL
cp $parafile $temppara

$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL2
cp Coverage/code/src/*/*.gc* coverage/testDL2
fi


#afterslip model with a single snapshot:
ln1="OutputForecastFile=output_testcases/testDL3"
ln2="Logfile=output_testcases/testDL3.log"
ln13="InversionStartDate=2004-09-28T17:00:24Z"	#so it's comparable to testDL5

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "13s+.*+$ln13+" | sed "s*slipmodelslist0*slipmodelslist1*" > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels3.dat" >> temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL3
cp Coverage/code/src/*/*.gc* coverage/testDL3
fi

#2 afterslips, with same geometry (single subfault)
ln1="OutputForecastFile=output_testcases/testDL4"
ln2="Logfile=output_testcases/testDL4.log"
ln7="InputCatalogFile=input_testcases/testDL/catalog2.dat"
ln13="InversionStartDate=2004-09-28T17:00:24Z"
ln6="ForecastEndDate=2005-10-15T17:15:24Z"
ln4="IssueDate=2005-09-30T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist2*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+" | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels4.dat" >> temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL4
cp Coverage/code/src/*/*.gc* coverage/testDL4
fi

#2 afterslip models, first one is empty (should be same as testDL3)
ln1="OutputForecastFile=output_testcases/testDL5"
ln2="Logfile=output_testcases/testDL5.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist1*" | sed "13s+.*+$ln13+" > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels5.dat" >> temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL5
cp Coverage/code/src/*/*.gc* coverage/testDL5
fi

#2 afterslip models, with splines
ln1="OutputForecastFile=output_testcases/testDL6"
ln2="Logfile=output_testcases/testDL6.log"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist6*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+"  | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels6.dat" >> temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL6
cp Coverage/code/src/*/*.gc* coverage/testDL6
fi

#1 afterslip models, corresponding to an empty earthquake. Should be the same as w/o the eqk (testDL8)
ln1="OutputForecastFile=output_testcases/testDL7"
ln2="Logfile=output_testcases/testDL7.log"
ln7="InputCatalogFile=input_testcases/testDL/catalog3.dat"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist3*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+"  | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels7.dat" >> temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL7
cp Coverage/code/src/*/*.gc* coverage/testDL7
fi

#1 afterslip models, corresponding to no earthquake. Should be the same as testDL7
ln1="OutputForecastFile=output_testcases/testDL8"
ln2="Logfile=output_testcases/testDL8.log"
ln7="InputCatalogFile=input_testcases/testDL/catalog4.dat"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  | sed "s*slipmodelslist0*slipmodelslist3*" | sed "13s+.*+$ln13+" | sed "7s+.*+$ln7+"  | sed "6s+.*+$ln6+" | sed "4s+.*+$ln4+" > temp_inputDL
echo "InputListAfterslipModels=input_testcases/testDL/aslipmodels7.dat" >> temp_inputDL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputDL
if [ $Build == "Coverage" ]
then
mkdir coverage/testDL8
cp Coverage/code/src/*/*.gc* coverage/testDL8
fi

