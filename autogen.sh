#!/bin/bash

temp=$(mktemp)

plot_out=mdc_$(./run_counter.sh).ps

cat <<EOF >>$temp

< garftrans
!add meta type PostScript file-name "$plot_out"
!open meta
!act meta

&CELL
solids
EOF
./generate_cell.py >>$temp
cat <<EOF >>$temp

nebem max-elem 5

opt nodebug
&FIELD
opt nodebug

// View the structure
area -1 -1 -1 1 1 1 ...
   nooutline ...
    light -20 -40 ...
    view -3*x+5*y-1*z=0 3d
Call plot_field_area
Call plot_end

# //  contour plot along wire
# area -0.4 -0.4 -0.4 0.4 0.4 0.4 ...
#     view  cut rot 90
# grid 25
# pl contour

EOF


## movie
#for i in $(seq -0.1 0.01 0.1); do
#cat <<EOF>>$temp
#//  contour plot
#area -1 -0.4 -0.4 1 0.4 0.4 ...
#    view x=$i cut rot 90
#grid 25
#pl contour
#
#EOF
#done





export LD_LIBRARY_PATH="./"
# cat $temp | ./garfield-9
# gv $plot_out
cat $temp