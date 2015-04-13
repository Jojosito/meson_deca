#!/usr/bin/env python
# plot_2d_csv.py
#
# NAME
#     plot_2d_csv.py - plot the two last parameters of the CmdSTAN output file
#
# SYNOPSIS
#     ./plot_2d_csv.py OUTPUT_CSV OUTPUT_PDF
#     ./plot_2d_csv.py
#
# DESCRIPTION
#     plot_2d_csv.py makes 1d histogram plots of the last and second to last
#     parameters stored in the CmdSTAN file OUTPUT_CSV (default : 'output.csv') 
#     and saves it to the file OUTPUT_PDF (default : 'output.pdf'). It also
#     saves the joint 2d histogram of the two last parameters.

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
print('Data shape:'.format(data.shape))

# Extract the last two variables from the data
theta_1_re = [x[-6] for x in data]
theta_1_im = [x[-5] for x in data]
theta_2_re = [x[-4] for x in data]
theta_2_im = [x[-3] for x in data]
theta_3_re = [x[-2] for x in data]
theta_3_im = [x[-1] for x in data]

print('Plotting the histograms...')
# Plot the data
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
nbins = 100

H1, xedges1, yedges1 = np.histogram2d(theta_1_re, theta_1_im, bins=nbins)
ax1.pcolor(yedges1, xedges1, H1)
H2, xedges2, yedges2 = np.histogram2d(theta_2_re, theta_2_im, bins=nbins)
ax2.pcolor(yedges2, xedges2, H2)
H3, xedges3, yedges3 = np.histogram2d(theta_3_re, theta_3_im, bins=nbins)
ax3.pcolor(yedges3, xedges3, H3)
ax4.hist(theta_2_re, bins=nbins)

plt.savefig(f_out_name)
plt.savefig(f_out_name[:-4] + '.png')






































