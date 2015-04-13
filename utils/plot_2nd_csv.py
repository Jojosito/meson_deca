#!/usr/bin/env python
# plot_csv.py
#
# NAME
#     plot_csv.py - plot the last parameter of the CmdSTAN output file
#
# SYNOPSIS
#     ./plot_csv.py OUTPUT_CSV OUTPUT_PDF
#     ./plot_csv.py
#
# DESCRIPTION
#     plot_csv.py makes a 1d histogram plot of the last parameter stored in the
#     CmdSTAN file OUTPUT_CSV (default : 'output.csv') and saves it to the file
#     OUTPUT_PDF (default : 'output.pdf').

import numpy as np
import matplotlib.pyplot as plt
import sys



if len(sys.argv) == 3:
    f_in_name = sys.argv[1]
    f_out_name = sys.argv[2]

elif len(sys.argv) == 1:
    f_in_name = 'output.csv'
    f_out_name = 'output.pdf'

else:
    sys.exit('Expected 0 or 2 arguments, got {0} instead. Aborting...'.format(len(sys.argv)))

#
print('plot_csv.py: Reading {0}...'.format(f_in_name))

f_in = open(f_in_name)

# Pack all numbers into 'data'
data = []
print('Importing data...')

for line in f_in:
    # Ignore lines containing STAN commands and STAN info...
    if line[0] in ['#', 'l', '\n', ' ']:
        pass
    else:
         # a is the list containing data of the line in the string form
         a = line.split(",")
         a = [float(x) for x in a if x != '']
         data.append(a)

data = np.asarray(data)
print('Data shape: {0}.'.format(data.shape))

# Extract the second to last variable from the data
m12 = [x[-2] for x in data]

print('Plotting the histograms...')
# Plot the data
f, ax1 = plt.subplots(1, 1)
nbins = 100
ax1.hist(m12, bins=nbins)
#ax2.plot(m12, '.')
plt.savefig(f_out_name)

print("plot_csv.py: Done. Plot saved in {0}.".format(f_out_name))

