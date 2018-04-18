#!/bin/bash




temp=$(mktemp)

plot_out=mdc_$(./run_counter.sh).ps
mkdir tracks

cell_generator=./MDC_cells/gen_MDC_II_40deg_1700V_3x3x3.py
gas_def=./gasses/ar70_co230.txt
# gas_def=./gasses/ar86_co214.txt

sense_angle_deg=40
cathode_pitch=0.2

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
drift_UV_charge=true

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
sin_angle=$(echo "s($sense_angle_deg/360*2*4*a(1))"| bc -l)
cos_angle=$(echo "c($sense_angle_deg/360*2*4*a(1))"| bc -l)
tan_angle=$(echo "$sin_angle/$cos_angle" | bc -l)
plane_x_coef=$(echo "1/$tan_angle" | bc -l)
cathode_repetition_halflength=$(echo "$cathode_pitch/$cos_angle/2" | bc -l)

for i in $(seq -$cathode_repetition_halflength 0.02 $cathode_repetition_halflength); do
plane_sum=$(echo "$i / $sin_angle" | bc -l)
area_x0=$(echo "-1   + $i*$cos_angle" | bc -l)
area_x1=$(echo " 1   + $i*$cos_angle" | bc -l)
area_y0=$(echo "-0.4 + $i*$sin_angle" | bc -l)
area_y1=$(echo " 0.4 + $i*$sin_angle" | bc -l)
area_z0=-0.4
area_z1=0.4
cat <<EOF>>$temp
//  contour plot
area $area_x0 $area_y0 $area_z0 $area_x1 $area_y1 $area_z1 ...
    view y+$plane_x_coef*x=$plane_sum cut rot 90
grid 25
pl surf
**pl cont

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
***track 0.1 -0.001 0.1  0.1 0.001 0.1 lines 230
track 0.1 -0.001 0.1  0.1 0.001 0.1 lines 20
track -0.1 -0.001 0.1  -0.1 0.001 0.1 lines 20
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
n=1

for i in $(seq 1 $n); do
cat <<EOF >>$temp

area -0.40 -0.40 -0.40 0.40 0.40 0.40 view x=0 3d
track $x $y -$halflength  $x $y $halflength HEED electron energy  0.70 GeV
INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
DRIFT TRACK MONTE-CARLO-DRIFT LINE-PRINT

EOF
done
fi

###  simulate UV charge ####
if [ $drift_UV_charge == "true" ]; then 
root_drift_times=true
point_list="./MDC_II_sim_Laser.data/point_list.txt"

export x_offset=0 # µm
export y_offset=0
export z_offset=0

export charges=230
export events=100
export sigma_y=0.001 # mm
export sigma_x=0.001 # mm
export sigma_z=0.20451 # mm -> rayleigh length of non-gaussian ionization density in z
export rayleigh=0.20451 # mm -> rayleigh length of non-gaussian ionization density in z
export n_photons=3 # 3 photon absorbtion
# export charges=230 # ca Fe55 charge deposition
export outfile="./tracks/input_tracks.txt"


while read line <&3; do

  list_x=$(echo $line | cut -f 1 -d " ") # in um
  list_y=$(echo $line | cut -f 2 -d " ") # in um
  list_z=$(echo $line | cut -f 3 -d " ") # in um
  
  export x=$(echo "($list_x  + $x_offset)/1000" | bc -l) # in mm
  export y=$(echo "($list_y  + $y_offset)/1000" | bc -l) # in mm
  export z=$(echo "($list_z  + $z_offset)/1000" | bc -l) # in mm
    
  echo "point:" $line
  rm ./tracks/input_tracks.txt
  root -l 'track_generators/gen_UV_charge.C()' -q ## takes input arguments from environmet variables, look at all "export" statements above
  cat ./tracks/input_tracks.txt >> $temp
    
    
done 3<$point_list



fi



######## finally execute garfield ##########



export LD_LIBRARY_PATH="./"
./garfield-9 < $temp | tee garfield_stdout.txt
rm $temp # clean up garfield input data
ps2pdf $plot_out
convert $plot_out -alpha off -delay 400 $plot_out.gif

if [ $root_drift_times ]; then
  rm ./tracks/*
  csplit -f './tracks/' -b '%06d' garfield_stdout.txt "/^  new UV event/" '{*}'
  rm track_drift_line_data.txt
  here=$(pwd)
  cd tracks
  for i in *; do
    # get laser position:
    laser_position=$(grep "  laser position:" $i| sed "s/  laser position://")
    perl -pi -e "s/^/$i  /g" $i
#     grep -P "Hit S" $i |  perl -pi -e "s/unavailable.*//g" >> track_drift_line_data.txt
    grep -P "Hit [ST]" $i |  perl -p -e "s/unavailable unavailable unavailable  Hit [ST] solid//g"| perl -p -e "s/^/$laser_position /g" >> track_drift_line_data.txt
  done
  cd $here
  root -l -q 'ascii_to_ttree_laser.C("tracks/track_drift_line_data.txt")' 
fi

if [ $animation == "true" ]; then
  rm $plot_out
  xdg-open $plot_out.gif
else 
  rm $plot_out
  gv $( echo $plot_out | sed s/.ps/.pdf/ )
fi

