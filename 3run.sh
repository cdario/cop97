#!/bin/bash

for i in $(seq 1 5) # 
do
    (./cop97) & 
    if(( $i % 5==0)); then wait; fi # 5 concurrent subshells at most
done
wait
