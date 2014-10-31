#! /bin/bash

for i in $(ls input_testcases --color=never | grep test); do input_testcases/$i/run.sh; done
