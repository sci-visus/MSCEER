#!/bin/bash

echo "to run: "
echo "       ./runIsoStudy.sh 1:filename(raw_blur_volume) 2:scaleFactor(optional, default 5) 3:isoValue(optional default 14000.0)  4:interactive(default 0) 5:X (optional, default extracted from filename) 6:Y(optional, ...) 7:Z(optional, ...) "
echo " "
echo "to run iso study manually:"
echo "./script  dir_to/vol_blur_file  X  Y  Z  scaleFactor isovalue interactive"
echo ""
now=$(date +"%T")
iso_log="isoStudy_at_$now.log"

echo "Created log file with terminal output: $iso_log"
echo "located in data directory"

script="/home/sam/Documents/PhD/Research/GradIntegrator/build/intensityFieldContourView/intensityFieldContourView"

cwd=$(pwd)
data_dir="/home/sam/Documents/PhD/Research/GradIntegrator/build/intensityFieldContourView"

iso_log="$data_dir/$iso_log"


echo "Log For Iso Study ran at time $now" >> $iso_log

vol_blur_file=$1
#"10D1_626_1024_1024.raw_blur.raw"

# Get dimensions from filename                                                  
X=$(echo "$vol_blur_file" | cut -d '.' -f 1 | cut -d '_' -f 2)
Y=$(echo "$vol_blur_file" | cut -d '.' -f 1 | cut -d '_' -f 3)
Z=$(echo "$vol_blur_file" | cut -d '.' -f 1 | cut -d '_' -f 4)
# or used passed param
X=${5:-"$X"}
Y=${6:-"$Y"}
Z=${7:-"$Z"}

echo "X: $X, Y: $Y, Z: $Z" > $iso_log

interactive=0
interactive=${4:-"$interactive"}

isovalue=0
isovalue=${3:-"$isovalue"}

scaleFactor=5
scaleFactor=${2:-"$scaleFactor"}

$script "$data_dir/$vol_blur_file"  $X $Y $Z  $scaleFactor $isovalue $interactive > $iso_log
