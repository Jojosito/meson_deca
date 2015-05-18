#!/bin/bash

# initialize_model.sh
#   This script installs a meson_deca model into a specified folder
#
# USAGE
#   initialize_model.sh FOLDER_NAME
#
# DESCRIPTION
#   Creates a folder model/FOLDER_NAME; copies STAN files
#   necessary for data generation and parameter fitting into this 
#   folder. Also, copies the python module 'model.so' that describes
#   A_cv, num_resonances, num_variables.


###### FUNCTIONS
function cdmeson_deca
{
  while [[ $PWD != '/' && ${PWD##*/} != 'meson_deca' ]]; do cd ..; done
}

###### MAIN
# Define model directory and meson_deca directory
MODEL_FOLDER=$1
cdmeson_deca
MESON_DECA=$(pwd)


cd $MESON_DECA
# Copy STAN files to the model directory
mkdir -p models/$MODEL_FOLDER
mkdir -p models/$MODEL_FOLDER/backup
cp $MESON_DECA/lib/stan_lib/* $MESON_DECA/models/$MODEL_FOLDER
cp $MESON_DECA/lib/c_lib/model.hpp $MESON_DECA/models/$MODEL_FOLDER/backup/model.hpp

# Build python library containing model description
echo "#### MESON_DECA: Building python modules..."
cd "$MESON_DECA/lib/c_lib/py_wrapper"
# Delete old modules
rm -r build/*
python setup.py build

# Copy it to the model directory
cp "$MESON_DECA"/lib/c_lib/py_wrapper/build/lib*/model.so "$MESON_DECA/models/$MODEL_FOLDER"

# DEPRECATED. Explicit is better than implicit - let the user calculate the
# integrals.
# Calculate the normalization integral for the model
# (We pass integration boundaries as arguments to utils script here,
# but maybe it would be better to store them in model.so? (i.e. in model.hpp))
#echo "#### MESON_DECA: Calculating normalization integrals..."
#cd "$MESON_DECA/models/$MODEL_FOLDER"
#"$MESON_DECA"/utils/calculate_normalization_integral.py -1 1 -1 1 #0 3 0 3
#"$MESON_DECA"/utils/calculate_gendata_normalization_integral.py -2 2 0 0 -2 2 0 0 -2 2 0 0

# Add calculated I_gen integral to the existing data file
#cat STAN_data_generator_I_gen.data.R >> STAN_data_generator.data.R
#rm STAN_data_generator_I_gen.data.R
echo "Done."
