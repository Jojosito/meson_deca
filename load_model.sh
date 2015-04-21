#!/bin/bash

# load_model.sh
#
# Copy the model from backup/model.hpp to lib/c_lib.
# MUST be called from a model folder

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

# Copy the model file
cp backup/model.hpp $MDECA_DIR/lib/c_lib
# Rebuild the python module
$MDECA_DIR/wrap_python.sh

echo "MESON_DECA: A new model has been loaded."
echo "Call utils/calculate_normalization_integral.py,"
echo "if the variable bounds (Dalitz plot bounds in 2D case)"
echo "has been changed."
