!add meta type PostScript file-name "mdc.ps"
!open meta
!act meta

*****************CELL***********************
&CELL
solids
wire  centre 0 0.0 0.0 ...
      direction 1 0 0 ...
      radius 0.001 ... 
      half-length 0.5 ...
      conductor-2 ...
      label s ...
      voltage 0
wire  centre 0 0.3 0.0 ...
      direction 1 0 0 ...
      radius 0.001 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre 0 -0.3 0.0 ...
      direction 1 0 0 ...
      radius 0.001 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre 0.1 0 0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre -0.1 0 0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre 0.1 0 -0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre -0.1 0 -0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre 0.3 0 0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre -0.3 0 0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre 0.3 0 -0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700
wire  centre -0.3 0 -0.25 ...
      direction 0 1 0 ...
      radius 0.004 ... 
      half-length 0.5 ...
      conductor-1 ...
      voltage -1700

nebem max-elem 5

*opt nodebug

*****************FIELD**********************
&field
*track FROM 0. 1. TO 0. -1.
*grid 40 40 

// area -0.4 -0.4 -0.4 0.4 0.4 0.4 view y=0
// grid 50
// plot-field surface e
// 
// area -0.4 -0.4 -0.4 0.4 0.4 0.4 view x=0
// grid 50
// plot-field surface e

// area -0.4 -0.4 -0.4 0.4 0.4 0.4 view x=0 rot 90
// plot-field contour


*****************MAGNETIC*******************
&magnetic
components  0.00  0.00  0.00 T 

*****************GAS************************
&gas
gas-id "ar/co2 70/30"
global gas_file `ar_co270-30-293-760.dat`
temperature 293
pressure 760
Call inquire_file(gas_file,exist)
If exist Then
Say "Gas file exists, retrieving ..."
get {gas_file}
Else
Say "Gas file not found, generating ..."
magboltz argon 70. carbon-dioxide 30. mobility 1.000000
write dataset "ar_co270-30-293-760.dat" gasdata remark "magboltz-ar70-co230"
Endif
heed argon 70. carbon-dioxide 30.

*****************DRIFT**********************
&drift
SEL s
// area -0.40 -0.40 -0.40 0.40 0.40 0.40 view -3*x+5*y-1*z=0 3d 
// plot cont
// track 0 0.1 -1 0 0.1 1 lines 25
// DRIFT TRACK l-pr
// track -0.2  -0.2 -0.05 -0.2  -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// area -0.4 -0.4 -0.4 0.4 0.4 0.4 view x=0 rot 90
// plot hist diffusion vector vdx, vdy surf cont

// ****positron drift lines****
// area -0.40 -0.40 -0.40 0.40 0.40 0.40 view x=0.1 cut rot 90
// drift wires lines 50 l-pr // why does it not work



// ***output electron speed at coordinates
// speed 0.1 0.1 0.1 ELECTRON
area -0.4 -0.4 -0.4 0.4 0.4 0.4 view x=0 rot 90 3d
// track -0.1  -0.1 -0.05 -0.1  -0.1 0.25 lines 10
// DRIFT TRACK LINE-PRINT ION POSITIVE
// DRIFT TRACK TIME-GRAPH
// DRIFT TRACK LINE-PRINT ELECTRON NEGATIVE MONTE-CARLO-DRIFT
// DRIFT TRACK LINE-PRINT ION NEGATIVE
// DRIFT TRACK LINE-PRINT ELECTRON POSITIVE
// xt-plot
area -0.40 -0.40 -0.40 0.40 0.40 0.40 view -3*x+5*y-1*z=0 3d 

// for i from -10 to 10 do
// Global x=i/100
// track {x} 0 -0.25 {x} 0 0.25 lines 50
// // DRIFT TRACK TIME-GRAPH NOLINE-PLOT
// DRIFT TRACK 
// enddo

// for i from -10 to 10 do
// Global y=i/100
// track 0 {y} -0.25 0 {y} 0.25 lines 50
// DRIFT TRACK TIME-GRAPH NOLINE-PLOT
// // DRIFT TRACK 
// enddo

track 0 -0.15 -0.25  0 -0.15 0.25 lines 500
INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
DRIFT TRACK TIME-GRAPH MONTE-CARLO-DRIFT LINE-PRINT


// // v curve analog to COSY beamtime
// track 0 -0.3 0  0 0.3 0 lines 300
// INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
// DRIFT TRACK TIME-GRAPH MONTE-CARLO-DRIFT
// 
// v curve analog to COSY beamtime
// track 0.1 -0.3 0  0.1 0.3 0 lines 500
// INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
// DRIFT TRACK TIME-GRAPH MONTE-CARLO-DRIFT LINE-PRINT

// 'Laser ionization tracklet'
// track 0.1 -0.001 0.1  0.1 0.001 0.1 lines 300
// INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500
// DRIFT TRACK TIME-GRAPH MONTE-CARLO-DRIFT LINE-PRINT


// * produces something
// arrival electron 5 last dataset "arrival/electron.5" thresh 0.8

// plot-field contour

// plot contour time range 0.1 0.3
// track -0.15 -0.2 -0.05 -0.15 -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// track -0.1  -0.2 -0.05 -0.1  -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// track -0.05 -0.2 -0.05 -0.05 -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// track  0.0  -0.2 -0.05  0.00 -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// track  0.05 -0.2 -0.05  0.05 -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// track  0.1  -0.2 -0.05  0.1  -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// track  0.15 -0.2 -0.05  0.15 -0.2 0.25 lines 12
// DRIFT TRACK l-pr
// track  0.2  -0.2 -0.05  0.2  -0.2 0.25 lines 12
// DRIFT TRACK l-pr

// *****************SIGNAL*********************
// &signal
// *Track starting from left side
// area -0.40 -0.40 -0.40 0.40 0.40 0.40 view y=0 rot 180 3d
// SELECT s
// AVALANCHE fixed 44130
// INTEGRATION-PARAMETERS  MONTE-CARLO-COLLISIONS 500 M-C-DIST-INT 0.000500 MAX-STEP 0.001000 INT-ACC 1.e-10 PROJECTED TRAP-RADIUS 1.010000
// OPTIONS NOCLUSTER-PLOT NOCLUSTER-PRINT
// ****************Loop over distance**********
// For j From 1 To 1 Do
// *-------------------------------------------
// Global r= 0.01*j
// Global alphaDEG=0.0
// Global alphaRAD=pi/180.*alphaDEG
// Global dx=0.1
// Global dy=dx*tan(pi/2.+alphaRAD)
// Global x1=r*cos(alphaRAD)
// Global y1=r*sin(alphaRAD)
// Global x=( 0.25-(y1-x1*tan(pi/2.+alphaRAD)))/tan(pi/2.+alphaRAD)
// Global y=  0.25
// TRACK FROM {x} {y} DIRECTION {dx} {dy}  RANGE 7.000000 HEED electron energy  0.70 GeV
// *-------------------------------------------
// TIME-WINDOW 0.0 0.001 1000
// ****************Loop over signal************
// For i From 1 To 1  Do
// SIGNAL AVALANCHE NOELECTRON-PULSE ION-TAIL RUNGE-KUTTA-DRIFT-LINES
// PLOT-SIGNALS
// *WRITE-SIGNALS DATASET "garfield_mdc-0_angle-0.txt" REMARK "p {j} n {i}" UNITS NANO-SECOND
// enddo
// ****************End Loop over signal*********
// enddo
// ****************End Loop over distance*******
// 
// 