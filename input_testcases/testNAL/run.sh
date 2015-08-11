!/bin/bash

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
basefile="input_testcases/testNAL/input.txt"
parafile="input_testcases/parameters.txt"
temppara="input_testcases/testNAL/temp_par.txt"

#run with all times before mainshock:
ln1="OutputForecastFile=output_testcases/testNAL1"
ln2="Logfile=output_testcases/testNAL1.log"
ln4="ForecastStartDate=2004-09-26T17:15:24Z"
ln5="ForecastEndDate=2004-09-27T17:15:24Z"
ln6="IssueDate=2004-09-25T17:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputNAL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputNAL
if [ $Build == "Coverage" ]
then
mkdir coverage/testNAL1
cp Coverage/code/src/*/*.gc* coverage/testNAL1
fi


#run with all times before mainshock (forecast, inversion overlap):
ln1="OutputForecastFile=output_testcases/testNAL2"
ln2="Logfile=output_testcases/testNAL2.log"
ln4="ForecastStartDate=2004-09-24T17:15:24Z"
ln5="ForecastEndDate=2004-09-25T17:15:24Z"
ln6="IssueDate=2004-09-25T17:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputNAL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputNAL
if [ $Build == "Coverage" ]
then
mkdir coverage/testNAL2
cp Coverage/code/src/*/*.gc* coverage/testNAL2
fi


#run with all times after mainshock:
ln1="OutputForecastFile=output_testcases/testNAL3"
ln2="Logfile=output_testcases/testNAL3.log"
ln4="ForecastStartDate=2004-09-30T17:15:24Z"
ln5="ForecastEndDate=2004-10-03T17:15:24Z"
ln6="IssueDate=2004-09-30T17:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputNAL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputNAL
if [ $Build == "Coverage" ]
then
mkdir coverage/testNAL3
cp Coverage/code/src/*/*.gc* coverage/testNAL3
fi


#run with all times after mainshock (earlier IssueTime):
ln1="OutputForecastFile=output_testcases/testNAL4"
ln2="Logfile=output_testcases/testNAL4.log"
ln4="ForecastStartDate=2004-09-30T17:15:24Z"
ln5="ForecastEndDate=2004-10-03T17:15:24Z"
ln6="IssueDate=2004-09-30T07:15:24Z"
ln11="InversionStartDate=2004-09-24T17:15:24Z"

sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+" | sed "4s+.*+$ln4+" | sed "5s+.*+$ln5+" | sed "6s+.*+$ln6+" | sed "11s+.*+$ln11+" > temp_inputNAL
cp  $parafile $temppara
$Build/CRS_3.0 temp_inputNAL
if [ $Build == "Coverage" ]
then
mkdir coverage/testNAL4
p Coverage/code/src/*/*.gc* coverage/testNAL4
fi

rm temp_inputNAL
rm $temppara