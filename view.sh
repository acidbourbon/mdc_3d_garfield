#!/bin/bash

generator=$1

temp=$(mktemp)

plot_out=mdc_$(./run_counter.sh).ps



### define plot output ###
cat <<EOF >>$temp
!add meta type PostScript file-name "$plot_out"
!open meta
!act meta

EOF

### generate cell
./MDC_cells/$generator >>$temp

cat <<EOF >>$temp
&FIELD
opt nodebug

// View the structure
area -1.2 -1.2 -1.2 1.2 1.2 1.2 ...
   nooutline ...
    light -20 -40 ...
    view -3*x+5*y-1*z=0 3d
Call plot_field_area
Call plot_end

//  surface plot along wire
area -0.4 -0.4 -0.4 0.4 0.4 0.4 ...
    view y=0 cut rot 90
grid 25
pl surface e

//  surface plot cross section
area -0.4 -0.4 -0.4 0.4 0.4 0.4 ...
    view x=0 cut rot 90
grid 25
pl surface e

//  surface plot cross section
area -0.6 -0.6 -0.6 0.6 0.6 0.6 ...
    view x=0 cut rot 90
grid 50
pl surface e

//  plot cross section
area -0.6 -0.6 -0.6 0.6 0.6 0.6 ...
    view x=0 cut rot 90
Call plot_field_area
Call plot_end

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
./garfield-9 < $temp
ps2pdf $plot_out
gv $( echo $plot_out | sed s/.ps/.pdf/ )