#!/bin/bash

# fit.sh NUM_SAMPLES
#
# Fit the model using NUM_SAMPLES and plot the result.

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

# Do fitting and plot the results
# (Sample 4 chains)
for i in {1..3}
do
  ./STAN_amplitude_fitting sample id=$i data file=STAN_amplitude_fitting.data.R output file=output$i.csv & #init=STAN_data_generator.data.R
done
