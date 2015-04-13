#!/bin/bash

# clear.sh
#   Builds
#
#    *    STAN_data_generator
#    *    STAN_amplitude_fitting
#
#   from the corresponding *.stan files in the current folder.


###### FUNCTIONS
function cdmeson_deca
{
  while [[ $PWD != '/' && ${PWD##*/} != 'meson_deca' ]]; do cd ..; done
}

###### MAIN
# Define model directory and meson_deca directory
MODEL_FOLDER=$(pwd)
cdmeson_deca
MESON_DECA=$(pwd)


cd $MESON_DECA
cd ..
make $MODEL_FOLDER/STAN_data_generator
make $MODEL_FOLDER/STAN_amplitude_fitting
cd $MODEL_FOLDER
