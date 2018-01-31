#!/bin/bash

old_ctr=$(cat run_counter.txt)
printf "%03d" $( echo "$old_ctr + 1" | bc ) | tee run_counter.txt
