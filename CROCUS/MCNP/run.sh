#!/bin/bash

# This script runs the input file on Linux

# Absolute location of the MCNP executable file
# (Change this path to the actual location of your MCNP executable)
EXE="/home/camden/School/My_Programs/MCNP-thermal-otf/MCNP_6.3_SOURCE/mcnp-src/mcnp-6.3.0-Source/mcnp6/build/mcnp6"

# Absolute location of the cross-section directory
# (Change this path to the actual location of your MCNP data directory)
# export DATAPATH="/home/camden/School/Data/MCNP"

# Number of threads used for parallel execution
NUM_THREADS=16

# Input file
INP="model.mcnp"

# OTF file
# OTF="/home/camden/School/Data/MCNP/MCNP_DATA/ENDF80SABOTF/tsl-HinH2O_OTF.h5"
# OTF="/home/camden/School/Data/MCNP/MCNP_DATA/ENDF80SABOTF/tsl-crystalline-graphite_OTF.h5"
# OTF="/home/camden/School/Data/MCNP/MCNP_DATA/ENDF80SABOTF/tsl-HinYH2_OTF.h5"

# Run simulation
# $EXE inp=$INP otf $OTF tasks $NUM_THREADS
$EXE inp=$INP tasks $NUM_THREADS
