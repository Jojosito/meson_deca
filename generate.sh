#!/bin/bash

# generate.sh NUM_SAMPLES
#
# Generate the data using STAN_data_generator executable, with
# NUM_SAMPLES events. The samples are converted to a '.root' file 
# and to model-dependent '.data.R' file for future analysis.
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

# Generate and plot data
./STAN_data_generator sample num_samples=$1 data file=STAN_data_generator.data.R output file=generated_data.csv
$MDECA_DIR/utils/plot_2d_csv.py generated_data.csv generated_data.pdf


# Create tree
$MDECA_DIR/utils/csv_to_root.py
$MDECA_DIR/utils/data_analysis__root_to_dataR.py
