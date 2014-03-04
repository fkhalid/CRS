#!/bin/bash

log="test3.log"

Release/CRS_2.0 clean
for i in $(seq 0 20); do Release/CRS_2.0 input/input_files/test3/input$i.dat; cp output/oldLL/LL.dat output/$i.LL; done

#Cumulative version:
Release/CRS_2.0 clean
Release/CRS_2.0 input/input_files/test3/inputCumu.dat

#This uses more revent available slip model (from 2012), just for curiousity:
Release/CRS_2.0 clean
Release/CRS_2.0 input/input_files/test3/inputCumu2.dat


cd output/test3
if [ -f $log ]
then
rm $log
fi
for i in $(seq 0 20); do echo $i"days:    " $(grep hash "$i"days.log)>> $log; done

n=`wc -l 20days_foremap.dat | awk '{print $1}'`
g5=`paste 20days_foremap.dat 20daysCumu_foremap.dat | awk '{print $9,"\t", $18,"\t", ($9-$18)/$9}' | awk -v toll=0.05 'sqrt($3*$3)>toll{print $0}' | wc -l`
g1=`paste 20days_foremap.dat 20daysCumu_foremap.dat | awk '{print $9,"\t", $18,"\t", ($9-$18)/$9}' | awk -v toll=0.01 'sqrt($3*$3)>toll{print $0}' | wc -l`
g10=`paste 20days_foremap.dat 20daysCumu_foremap.dat | awk '{print $9,"\t", $18,"\t", ($9-$18)/$9}' | awk -v toll=0.1 'sqrt($3*$3)>toll{print $0}' | wc -l`

echo "">>$log
echo "Comparing files 20days_foremap.dat 20daysCumu_foremap.dat:" >> %log
echo $g1" out of "$n" gridpoints different by more than 1%" >> $log
echo $g5" out of "$n" gridpoints different by more than 5%" >> $log
echo $g10" out of "$n" gridpoints different by more than 10%" >> $log

cd ../..
