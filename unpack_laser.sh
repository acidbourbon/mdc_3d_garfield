#!/bin/bash
infile=$1
root -l "ascii_to_ttree_laser.C(\"$infile\")" $2 $3 $4 $5 $6
