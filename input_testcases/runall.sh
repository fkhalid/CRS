#! /bin/bash

if [ $# -eq 0 ]
then
 Build="Release"
else
 Build=$1
fi

if [ $Build == "Coverage" ]
#if [ 0 -eq 1 ]
then
   if [ -e coverage ]
   then
	echo -n "Overwrite folder coverage [y/n]?"
	read answer

	if [ $answer == "y" ]
	then
	echo "Will overwrite folder."
	rm -rf coverage
	elif [ $answer == "n" ]
	then
	echo -n "How should old 'coverage' be renamed?"
	read newname
	echo "Will move coverage to $newname"
	mv coverage $newname
	else
	echo "Invalid answer."
	fi
   fi

   mkdir coverage
fi

#Run code:
for i in $(ls input_testcases --color=never | grep test); do input_testcases/$i/run.sh $Build; done
#for i in $(ls input_testcases --color=never | grep 'testF'); do input_testcases/$i/run.sh $Build; done

#---------------------------------#
#Put together coverage information
#---------------------------------#
if [ $Build == "Coverage" ]
then
#copy coverage folder to specific one for this run:
covdir=` echo "coverage_testcases/"$(date +%d_%m_%Y_%H%M%S)`
mv coverage $covdir
cp -r code $covdir

#run gcov on each folder to get partial coverage:
cd $covdir

for f in $(ls --color=never| grep test);
   do

   cd $f
   for i in $(ls ../code/src/*/*.c)
	do 
	gcov-4.4 $i -o .
	done
   #produce individual html file:
   gcovr --html -g -o gcov.html -k
   cd ..
done

#merge gcov files to get a cumulative files:
mkdir cumulative
for i in $(ls --color=never $(ls --color=never | grep test | head -1) | grep gcov)
   do
   gcov_merge.sh $(ls */$i | grep -v cumulative) > cumulative/$i
   gcov_simplify.sh cumulative/$i
   done

cd cumulative

#produce cumulative html file:
gcovr --html -g -o gcov.html -k
cd ../..

fi

