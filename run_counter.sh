#!/bin/bash

if [ -e "run_counter.txt" ]; then
old_ctr=$(cat run_counter.txt)
printf "%03d" $( echo "$old_ctr + 1" | bc ) | tee run_counter.txt
else
printf "%03d" 0 | tee run_counter.txt
fi