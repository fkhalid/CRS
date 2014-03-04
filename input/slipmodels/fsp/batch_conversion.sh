#!/bin/bash

for i in $(ls | grep fsp)
do
n=`echo $i | sed 's-.fsp--'`
echo $i $n;
~/Code/slipmodels_conversions/bin/fsp2pscmp $i ../inp/$n
/home/des/camcat/Code/Scripts/slipmodel_sep.sh ../inp/$n
done

