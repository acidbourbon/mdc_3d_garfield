#!/bin/bash



mkdir -p graphics_out
mkdir -p tracks



temp=$(mktemp)

plot_out=graphics_out/mdc_$(./run_counter.sh).ps

  rm ./tracks/*

cell_generator=./MDC_cells/gen_MDC_I_1600V_3x3x3.py
# gas_def=./gasses/ar70_co230.txt
gas_def=./gasses/ar86_co214.txt


# animation=true
# root_drift_times=true


##### the switches #####

# view_cell=true
# plot_field=true
# field_animation=true
# view_inner_cube=true
# drift_point_charge=true
# drift_MIPS_vcurve=true
# drift_MIPS_track=true
drift_cosmics=true

cat garfinit.txt >> $temp



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
if [ $view_cell == "true" ]; then cat <<EOF >>$temp

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
if [ $view_inner_cube == "true" ]; then cat <<EOF >>$temp
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
if [ $plot_field == "true" ]; then cat <<EOF >>$temp

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



### no magnetic field
cat <<EOF>>$temp
*****************MAGNETIC*******************
&magnetic
components  0.00  0.00  0.00 T 

EOF


## gas section
cat $gas_def >>$temp

### drift calculations
cat <<EOF>>$temp
******************DRIFT*********************
&drift
SEL s

EOF


###  simulate drift of point charge ####
if [ $drift_point_charge == "true" ]; then 
root_drift_times=true

cat <<EOF >>$temp

**** something like Fe55 point charge 
area -0.40 -0.40 -0.40 0.40 0.40 0.40
track 0.1 -0.001 0.1  0.1 0.001 0.1 lines 230
INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
DRIFT TRACK TIME-GRAPH MONTE-CARLO-DRIFT LINE-PRINT

EOF
fi

###  simulate MIPS v curve (as recorded at COSY) ####
if [ $drift_MIPS_vcurve == "true" ]; then 
root_drift_times=true

for y in $(seq -0.3 0.05 0.3); do

cat <<EOF >>$temp

area -0.40 -0.40 -0.40 0.40 0.40 0.40
track 0 $y -0.25  0 $y 0.25 HEED electron energy  0.70 GeV
INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
DRIFT TRACK MONTE-CARLO-DRIFT LINE-PRINT NOLINE-PLOT

EOF

done

fi

###  simulate single MIPS track) ####
if [ $drift_MIPS_track == "true" ]; then 
root_drift_times=true

x=0.1
y=0.1
halflength=0.4
n=3

for i in $(seq 1 $n); do
cat <<EOF >>$temp

area -0.40 -0.40 -0.40 0.40 0.40 0.40 view x=0 3d
track $x $y -$halflength  $x $y $halflength HEED electron energy  0.70 GeV
INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
DRIFT TRACK MONTE-CARLO-DRIFT LINE-PRINT

EOF
done
fi

###  simulate cosmics coming from above ####
if [ $drift_cosmics == "true" ]; then 
root_drift_times=true


cat <<EOF >>$temp

# area -0.40 -0.40 -0.40 0.40 0.40 0.40 view x=0 3d
# INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500

area -0.4 -0.4 -0.4 0.4 0.4 0.4 ...
   nooutline ...
    light -20 -40 ...
    view -3*x+5*y-1*z=0 3d
EOF

# void gen_cosmic_tracks(
#   TString outfile_str     = "input_tracks.txt",
#   TString number_str     = "20",
#   TString z_length_str     = "0.6",
#   TString displacement_x_str = "0.2",
#   TString displacement_y_str = "1.0"
# )

#http://pdg.lbl.gov/2011/reviews/rpp2011-rev-cosmic-rays.pdf
root -l 'track_generators/gen_cosmic_tracks.C("./tracks/input_tracks.txt","100000","0.6","0.1","0.7")' -q
cat ./tracks/input_tracks.txt >> $temp

fi



######## finally execute garfield ##########



export LD_LIBRARY_PATH="./"
./garfield-9 < $temp | tee garfield_stdout.txt
# ps2pdf $plot_out
# convert $plot_out -alpha off -delay 400 $plot_out.gif

if [ $root_drift_times ]; then
  csplit -f './tracks/' -b '%05d' garfield_stdout.txt "/^1 Track drift line plot :/" '{*}'
  rm track_drift_line_data.txt
  here=$(pwd)
  cd tracks
  for i in *; do
    perl -pi -e "s/^/$i  /g" $i
    grep -P "Hit S" $i |  perl -pi -e "s/unavailable.*//g" >> track_drift_line_data.txt
  done
  cd $here
  xterm -e root -l 'ascii_to_ttree.C("tracks/track_drift_line_data.txt")' &
fi

if [ $animation == "true" ]; then
  rm $plot_out
  xdg-open $plot_out.gif
else 
#   rm $plot_out
#   gv $( echo $plot_out | sed s/.ps/.pdf/ )
fi

