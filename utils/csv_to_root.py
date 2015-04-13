#!/usr/bin/env python
# csv_to_root.py

# NAME
#    csv_to_root.py - make a ROOT tree from the 'output.csv' STAN file
#
# SYNOPSIS
#    ./csv_to_root.py [OPTIONS...]
#
# DESCRIPTION
#    csv_to_root.py parses the specified file (by default, 
#    'generated_data.csv'), ignoring all lines starting with '#', and writes 
#    the entries of the file in a tree.

import argparse
import numpy as np
#import matplotlib.pyplot as plt
import os
import ROOT


# Parse the arguments
parser = argparse.ArgumentParser(description='Script to convert STAN output.csv to a .ROOT file.')

parser.add_argument('--tree_name', 
                    default='t',
                    help="Name of the tree in the .ROOT file.")

parser.add_argument('--tree_title', 
                    default='t',
                    help="Title of the tree in the .ROOT file.")

parser.add_argument('f_in', 
                     default='generated_data.csv',
                     nargs='?',
                     type=argparse.FileType('r'))

parser.add_argument('f_out', 
                    default='generated_data.root',
                    nargs='?',
                    type=argparse.FileType('w')
)

parser.add_argument('-p', '--print_dalitz_plot',
                    default=0,
                    type=float)

args = parser.parse_args()


#
print('csv_to_root.py: Reading {0}...'.format(args.f_in.name))

# Open CSV file
#print('Importing data...')
f_in = args.f_in

# Open ROOT file
f_out = ROOT.TFile(args.f_out.name, "recreate")
t = ROOT.TTree(args.tree_name, args.tree_title)


for line in f_in:
    # Ignore lines containing STAN commands and STAN info...
    if line[0] in ['#', '\n', ' ']:
        pass

    # Parse the line containing names of the parameters...
    elif line[0] == 'l':
        param_names = line.split(",")
        param_names = [x.replace("\n", "") for x in param_names if x != '']
        # Parameter values will be stored in param[0], param[1], ... etc.
        param = [np.zeros(1, dtype=float) for i in range(len(param_names))]
        param = np.asarray(param)
        for i in range(len(param_names)):
            t.Branch(param_names[i], param[i], param_names[i] + '/D')

    # Parse the lines containing parameter values...
    else:
         a = line.split(",")
         a = [float(x) for x in a if x != '']
         param[:,0] = [x for x in a]
         t.Fill()

# Make a *.pdf drawing of the parameters 'm2_ab', 'm2_bc'
num_bins = str(200)
if args.print_dalitz_plot == 1:
   c = ROOT.TCanvas("y.1:y.2", "Dalitz plot")
   t.Draw("y.1:y.2>>hh(" + num_bins + ", 0, 3.1, " + num_bins + ", 0, 3.1)", "", "COLZ")
   c.Print(args.f_out.name[:-5] + '.pdf' )


f_out.Write()
f_out.Close()

f_in.close()


print("csv_to_root.py: Done. Data saved in {0}.".format(args.f_out.name))
