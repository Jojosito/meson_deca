#!/bin/bash

# fit.sh
#
# Prints the file output.csv. Should be run after 'fit.sh' is finished.

###### FUNCTIONS
function cdmeson_deca
{
  while [[ $PWD != '/' && ${PWD##*/} != 'meson_deca' ]]; do cd ..; done
}

###### MAIN
# Define locations of CmdStan and current folder
MODEL_DIR=$PWD
cdmeson_deca
MDECA_DIR=$PWD
cd $MODEL_DIR

# Merge chain outputs
grep lp__ output1.csv > output.csv
sed '/^[#l]/d' output?.csv >> output.csv

# Plot results
$MDECA_DIR/utils/plot_6d_csv.py