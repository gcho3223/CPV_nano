#! /bin/bash

### Configuration ###

version="SelLep"
period="UL2016PreVFP"
channel="MuMu"

### Running ###
echo "root -l -b -q Plotter_Kinematics.C'(\"$version\",\"$period\",\"$channel\")'"
eval "root -l -b -q Plotter_Kinematics.C'(\"$version\",\"$period\",\"$channel\")'"