#!/bin/bash


temp=$(mktemp)

plot_out=mdc_$(./run_counter.sh).ps


cell_generator=./MDC_cells/gen_MDC_II_1700V_3x3x3.py

# view_cell=true
# plot_field=true
# field_animation=true
view_inner_cube=true




### define plot output ###
cat <<EOF >>$temp
!add meta type PostScript file-name "$plot_out"
!open meta
!act meta

EOF

### generate cell
$cell_generator >>$temp

### field section ###
cat <<EOF >>$temp
&FIELD
*opt nodebug

EOF

### view cell plots ###
if [ $view_cell == "true" ]; then
cat <<EOF >>$temp

// View the structure
area -1.2 -1.2 -1.2 1.2 1.2 1.2 ...
   nooutline ...
    light -20 -40 ...
    view -3*x+5*y-1*z=0 3d
Call plot_field_area
Call plot_end

//  plot cross section
area -0.6 -0.6 -0.6 0.6 0.6 0.6 ...
    view x=0.1 cut rot 90
Call plot_field_area
Call plot_end

EOF
fi

### view inner cube ###
if [ $view_inner_cube == "true" ]; then
cat <<EOF >>$temp
// View innermost cube
area -0.4 -0.4 -0.4 0.4 0.4 0.4 ...
   nooutline ...
    light -20 -40 ...
    view -3*x+5*y-1*z=0 3d
Call plot_field_area
Call plot_end

EOF
fi

### field plots ###
if [ $plot_field == "true" ]; then
cat <<EOF >>$temp

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


EOF
fi

## field animation
if [ $field_animation == "true" ]; then
animation=true
for i in $(seq -0.1 0.02 0.1); do
cat <<EOF>>$temp
//  contour plot
area -1 -0.4 -0.4 1 0.4 0.4 ...
    view x=$i cut rot 90
grid 25
pl surf

EOF
done
fi



## no magnetic field
cat <<EOF>>$temp
*****************MAGNETIC*******************
&magnetic
components  0.00  0.00  0.00 T 

EOF





######## finally execute garfield ##########



export LD_LIBRARY_PATH="./"
./garfield-9 < $temp | tee garfield_stdout.txt
ps2pdf $plot_out
if [ $animation == "true" ]; then
convert $plot_out -alpha off -delay 400 $plot_out.gif
xdg-open $plot_out.gif
else 
gv $( echo $plot_out | sed s/.ps/.pdf/ )
fi

