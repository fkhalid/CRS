#!/bin/bash

grep gcc -rl Debug | xargs sed -i "s*gcc*gcc -g -pg -Wall*g"
grep "\-O3" -rl Debug | xargs sed -i "s*-O3*-O0*g"
