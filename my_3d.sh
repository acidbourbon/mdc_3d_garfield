#!/bin/bash
export LD_LIBRARY_PATH=./
./garfield-9 < garfield_mdc-0_angle-0_3D.C | tee garfield_stdout.txt


# tdlp_line=$(grep -n "Track drift line plot :" garfield_stdout.txt | cut -f 1 -d ":")
# endprg_line=$(grep -n "  ------ INPGET MESSAGE : EOF on standard input ; end of program execution." garfield_stdout.txt | cut -f 1 -d ":")
# 
# head -n +$[ $endprg_line -1 ] garfield_stdout.txt | tail -n +$[ $tdlp_line + 13 ] | grep -P "^   \d" > track_drift_line_data.txt


grep -P "Hit S" garfield_stdout.txt |  perl -pi -e "s/unavailable.*//g" > track_drift_line_data.txt

ps2pdf mdc.ps; gv mdc.pdf 

root -l 'ascii_to_ttree.C("track_drift_line_data.txt")'