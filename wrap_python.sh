#!/bin/bash

# wrap_python.sh
#   This script rebuilds lib/c_lib/model.cpp to a python module and copies
#   that module into the current folder


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

cd "$MESON_DECA/lib/c_lib"
rm -r build
python setup.py build
cp "$MESON_DECA"/lib/c_lib/build/lib*/model.so "$MODEL_FOLDER"
cd "$MODEL_FOLDER"
