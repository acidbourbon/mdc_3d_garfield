#!/bin/bash

xlist=$(seq 0 0)
ylist="$(seq 0 100 3000)"
zlist=$(seq 0 0)

for x in $xlist; do
for y in $ylist; do
for z in $zlist; do

echo "$x $y $z"

done; done; done
