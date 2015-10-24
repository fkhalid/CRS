#!/bin/bash

ls code/src/*/* --color=never | grep -v "okada.c\|okada.h" | xargs sed -i '1s/^/\n/' 		#add empty line 'cause next line can only add stuff after a certain line.
ls code/src/*/* --color=never | grep -v "okada.c\|okada.h" | xargs sed -i '1r licence_stuff/headers.txt'	#add header contents

