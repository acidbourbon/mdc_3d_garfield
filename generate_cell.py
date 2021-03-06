#!/usr/bin/python

import numpy as np

zpitch=0.25
xpitch=0.2 # cathode pitch
ypitch=0.3 # sense/field wire pitch

field_wire_rad = 0.004
cath_wire_rad  = 0.004
sense_wire_rad = 0.001

v_cath=-1700 # cathode voltage
v_field=-1700 # field wire voltage
v_sense=0 # field wire voltage

sim_vol_xsize=1 # along the sense wire
sim_vol_ysize=1 # along the cathode wires
sim_vol_zsize=1 # in layer stacking direction

cath_conductor  = "conductor-1"
field_conductor = "conductor-1"
sense_conductor = "conductor-2"

## generate sense wires
x=0
for y in np.arange(0,sim_vol_ysize/2.,ypitch*2) :
  for z in np.arange(0,sim_vol_zsize/2.,zpitch*2) :
    if ( (z/zpitch)%4 ):
      y = y+ypitch
    print "wire  centre {:} {:} {:} ...".format(x,y,z)
    print "      direction 1 0 0 ..."
    print "      radius {:} ... ".format(sense_wire_rad)
    print "      half-length {:} ...".format(sim_vol_xsize/2.)
    print "      "+sense_conductor+" ..."
    print "      voltage {:}".format(v_sense)
    if (y > 0):
      print "wire  centre {:} {:} {:} ...".format(x,-y,z)
      print "      direction 1 0 0 ..."
      print "      radius {:} ... ".format(sense_wire_rad)
      print "      half-length {:} ...".format(sim_vol_xsize/2.)
      print "      "+sense_conductor+" ..."
      print "      voltage {:}".format(v_sense)
    if (z > 0):
      print "wire  centre {:} {:} {:} ...".format(x,y,-z)
      print "      direction 1 0 0 ..."
      print "      radius {:} ... ".format(sense_wire_rad)
      print "      half-length {:} ...".format(sim_vol_xsize/2.)
      print "      "+sense_conductor+" ..."
      print "      voltage {:}".format(v_sense)
    if (y > 0 and z > 0):
      print "wire  centre {:} {:} {:} ...".format(x,-y,-z)
      print "      direction 1 0 0 ..."
      print "      radius {:} ... ".format(sense_wire_rad)
      print "      half-length {:} ...".format(sim_vol_xsize/2.)
      print "      "+sense_conductor+" ..."
      print "      voltage {:}".format(v_sense)
  
## generate field wires
x=0
for y in np.arange(ypitch,sim_vol_ysize/2.,ypitch*2) :
  for z in np.arange(0,sim_vol_zsize/2.,zpitch*2) :
    if ( (z/zpitch)%4 ):
      y = y+ypitch
    print "wire  centre {:} {:} {:} ...".format(x,y,z)
    print "      direction 1 0 0 ..."
    print "      radius {:} ... ".format(sense_wire_rad)
    print "      half-length {:} ...".format(sim_vol_xsize/2.)
    print "      "+field_conductor+" ..."
    print "      voltage {:}".format(v_field)
    if (y > 0):
      print "wire  centre {:} {:} {:} ...".format(x,-y,z)
      print "      direction 1 0 0 ..."
      print "      radius {:} ... ".format(sense_wire_rad)
      print "      half-length {:} ...".format(sim_vol_xsize/2.)
      print "      "+field_conductor+" ..."
      print "      voltage {:}".format(v_field)
    if (z > 0):
      print "wire  centre {:} {:} {:} ...".format(x,y,-z)
      print "      direction 1 0 0 ..."
      print "      radius {:} ... ".format(sense_wire_rad)
      print "      half-length {:} ...".format(sim_vol_xsize/2.)
      print "      "+field_conductor+" ..."
      print "      voltage {:}".format(v_field)
    if (y > 0 and z > 0):
      print "wire  centre {:} {:} {:} ...".format(x,-y,-z)
      print "      direction 1 0 0 ..."
      print "      radius {:} ... ".format(sense_wire_rad)
      print "      half-length {:} ...".format(sim_vol_xsize/2.)
      print "      "+field_conductor+" ..."
      print "      voltage {:}".format(v_field)
  
  
## generate cathodes
y=0
for x in np.arange(xpitch/2.,sim_vol_xsize/2.,xpitch) :
  for z in np.arange(zpitch,sim_vol_zsize/2.,zpitch*2) :
    print "wire  centre {:} {:} {:} ...".format(x,y,z)
    print "      direction 0 1 0 ..."
    print "      radius {:} ... ".format(cath_wire_rad)
    print "      half-length {:} ...".format(sim_vol_ysize/2.)
    print "      "+cath_conductor+" ..."
    print "      voltage {:}".format(v_cath)
    if (x > 0):
      print "wire  centre {:} {:} {:} ...".format(-x,y,z)
      print "      direction 0 1 0 ..."
      print "      radius {:} ... ".format(cath_wire_rad)
      print "      half-length {:} ...".format(sim_vol_ysize/2.)
      print "      "+cath_conductor+" ..."
      print "      voltage {:}".format(v_cath)
    if (z > 0):
      print "wire  centre {:} {:} {:} ...".format(x,y,-z)
      print "      direction 0 1 0 ..."
      print "      radius {:} ... ".format(cath_wire_rad)
      print "      half-length {:} ...".format(sim_vol_ysize/2.)
      print "      "+cath_conductor+" ..."
      print "      voltage {:}".format(v_cath)
    if (x > 0 and z > 0):
      print "wire  centre {:} {:} {:} ...".format(-x,y,-z)
      print "      direction 0 1 0 ..."
      print "      radius {:} ... ".format(cath_wire_rad)
      print "      half-length {:} ...".format(sim_vol_ysize/2.)
      print "      "+cath_conductor+" ..."
      print "      voltage {:}".format(v_cath)


