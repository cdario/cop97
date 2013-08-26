#!/bin/bash

echo "Before starting subshell"
( #run tasks in parallel, different cores if any 

    count=1
    while [ $count -le 10 ]
    do
	echo "$count"
	./cop97 data0.txt 2 0
	sleep 1
	(( count++ ))
    done
) &
(
    count=1
    while [ $count -le 10 ]
    do
	echo "$count"
	./cop97 data0.txt 2 0
	sleep 1
	(( count++ ))
    done
) &
echo "Finished"
