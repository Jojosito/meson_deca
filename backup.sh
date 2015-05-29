#!/bin/bash

# backup.sh
#
# Copy the STAN*.stan and STAN*.data.R as well as
# lib/c_lib/model.hpp files to the backup folder
#
# CAVEAT: run from the model folder.

###### FUNCTIONS
function cdmeson_deca
{
  while [[ $PWD != '/' && ${PWD##*/} != 'meson_deca' ]]; do cd ..; done
}

###### MAIN
# Define locations of necessary files
MODEL_DIR=$PWD
cdmeson_deca
MDECA_DIR=$PWD
cd $MODEL_DIR

mkdir -p backup
cp STAN*.stan backup
cp STAN_data_generator.data.R backup
cp $MDECA_DIR/lib/c_lib/model.hpp backup

