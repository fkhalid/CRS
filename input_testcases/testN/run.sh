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
basefile="input_testcases/testN/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testN/temp_par.txt"

#run with all times before mainshock:
ln1="OutputForecastFile=output_testcases/testN1"
ln2="Logfile=output_testcases/testN1.log"
ln4="ForecastStartDate=2004-09-26T17:15:24Z"
ln5="ForecastEndDate=2004-09-27T17:15:24Z"
ln6="IssueDate=2004-09-25T17:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputN
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputN
if [ $Build == "Coverage" ]
then
mkdir coverage/testN1
cp Coverage/code/src/*/*.gc* coverage/testN1
fi


#run with all times before mainshock (forecast, inversion overlap):
ln1="OutputForecastFile=output_testcases/testN2"
ln2="Logfile=output_testcases/testN2.log"
ln4="ForecastStartDate=2004-09-24T17:15:24Z"
ln5="ForecastEndDate=2004-09-25T17:15:24Z"
ln6="IssueDate=2004-09-25T17:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputN
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputN
if [ $Build == "Coverage" ]
then
mkdir coverage/testN2
cp Coverage/code/src/*/*.gc* coverage/testN2
fi


#run with all times after mainshock:
ln1="OutputForecastFile=output_testcases/testN3"
ln2="Logfile=output_testcases/testN3.log"
ln4="ForecastStartDate=2004-09-30T17:15:24Z"
ln5="ForecastEndDate=2004-10-03T17:15:24Z"
ln6="IssueDate=2004-09-30T17:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputN
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputN
if [ $Build == "Coverage" ]
then
mkdir coverage/testN3
cp Coverage/code/src/*/*.gc* coverage/testN3
fi


#run with all times after mainshock (earlier IssueTime):
ln1="OutputForecastFile=output_testcases/testN4"
ln2="Logfile=output_testcases/testN4.log"
ln4="ForecastStartDate=2004-09-30T17:15:24Z"
ln5="ForecastEndDate=2004-10-03T17:15:24Z"
ln6="IssueDate=2004-09-30T07:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputN
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputN
if [ $Build == "Coverage" ]
then
mkdir coverage/testN4
cp Coverage/code/src/*/*.gc* coverage/testN4
fi

rm temp_inputN
rm $temppara
