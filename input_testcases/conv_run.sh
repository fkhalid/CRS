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
basefile="input_testcases/convergence_input.txt"
parafile="input_testcases/convergence_para.txt"
temppara="conv_tmppara.txt"

for i in $(seq 5 5)

do
 d=`echo "10^"$i | bc`
 ln1="OutputForecastFile=output_testcases/convB$d"
 ln2="Logfile=output_testcases/convB$d.log"

 sed "1s+.*+$ln1+" $basefile | sed "2s+.*+$ln2+"  > conv_temp_input
 sed '59s+.*+'$d'+' $parafile > $temppara
 $Build/run_crs conv_temp_input

done

rm conv_temp_input
rm $temppara
