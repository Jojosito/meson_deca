#!/bin/bash

# fit.sh
#
# Prints the file output.csv. Should be run after 'fit.sh' is finished.


# Merge chain outputs
grep lp__ output1.csv > output.csv
sed '/^[#l]/d' output?.csv >> output.csv

# Plot results
$MDECA_DIR/utils/plot_6d_csv.py