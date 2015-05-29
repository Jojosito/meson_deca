#!/usr/bin/env python

# data_analysis__root_to_dataR.py
#
# NAME
#    root_to_camp_dataR.py - convert ROOT tree to (model-dep) PWA amplitudes.
#
# DESCRIPTION
#    The events for our model are generated as a real-valued vector
#    y[] (for example, a vector containing invariant mass
#    squares). However, when we do the sampling over our events, it is
#    convenient to work with an array A_r(y) containing
#    model-dependent PWA amplitudes of our events. With other words,
#    to each event y_i = (m2_ab(i), m2_bc(i), ...), we would like to
#    compute complex numbers A_1(y_i), A_2(y_i), ... A_R(y_i) and
#    store them in a file that STAN can read - with other words, as
#    an array A_r_data dumped in the file
#    STAN_amplitude_fitting.data.R (or in the f_out argument).
#
#    To sum up, this script: opens a .ROOT file containing a Tree
#    with branches y.1 ... y.R (R meaning the number of variables);
#    gets the model-dependent amplitude function A_r(y) from model.so;
#    applies it to all events in the tree and dumps the result to 
#    *.data.R file.
#
# USAGE
#    data_analysis__root_to_dataR.py f_in f_out
#
#    Takes the ROOT file f_in ('generated_data.root' by default),
#    saves results in f_out ('STAN_amplitude_fitting.data.R' by default).
#
# CAVEAT
#    MUST be called from the model folder containing the .root file
#    and where the output file will be saved.

import argparse
import numpy as np
import os
import sys

MODEL_FOLDER = os.getcwdu()
sys.path.insert(1, MODEL_FOLDER)
import model # Contains A_r - model-dependent PWA resonance function

# Next line is bad - works only for models in the folder 'models', not deeper
sys.path.insert(1, "../../lib/py_lib")
import convert # Translates A_r results to usable form
from stan_rdump import *


# Parse the arguments
parser = argparse.ArgumentParser(description='Script to convert a .ROOT file to STAN data.R file.')


parser.add_argument('f_in',
                     default='generated_data.root',
                     nargs='?',
                     type=argparse.FileType('r'))

parser.add_argument('f_out',
                    default='STAN_amplitude_fitting.data.R',
                    nargs='?',
                    type=argparse.FileType('w')
)


args = parser.parse_args()

# User must know what's happening!
print("data_analysis__root_to_dataR.py: Reading {0}.".format(args.f_in.name))

# Load the branches y from *.root into python cache
f_in = ROOT.TFile(MODEL_FOLDER + '/' + args.f_in.name)
t = f_in.Get("t")

# Declare the variables that will contain tree entries
# Caveat: model already tells us how much variables we need
y = [np.asarray(0, dtype=float) for i in range(model.num_variables())]

# Aliase the variables to the corresponding branches
for i in range(model.num_variables()):
    # Branch names inherited from STAN start from 1: y.1, y.2 etc.
    # If the tree is not generated by STAN, but contains actual data,
    # it is important to rename the trees to y.1, y.2 ... - or to
    # change the code in this for-loop.
    branch_name = "y.{0}".format(i+1)
    t.SetBranchAddress(branch_name, y[i])

# D_ is the number of events
D_ = t.GetEntries()
y_data_ = np.asarray([np.zeros(D_, dtype = float) for i in range(model.num_variables())])


### FILL DATA FROM TREE ###
for d in range(D_):
    t.GetEntry(d)
    y_data_[:,d] = y


# Evaluate A_cv_ at y_data_
A_cv_data_ = np.asarray([convert.MatrixForm(model.A_cv(model.num_variables(), y_data_[:,d].tolist())) for d in range(D_)])

# Evaluate A_v_background_abs2_data_ at y_data_
A_v_background_abs2_data_ = np.asarray([convert.VectorForm(model.A_v_background_abs2(model.num_background(), y_data_[:,d].tolist())) for d in range(D_)])

# Define the integrals for the normalization function
# Usually these integrals can be generated by calling
# utils/calculate_normalization_integral.py from the model
# folder.
execfile('normalization_integral.py')
# Convert complex array to 2d-float array
I_ = np.transpose(I_)
I_out_ = np.zeros([2, model.num_resonances(), model.num_resonances()])
I_out_[0,:,:] = I_.real
I_out_[1,:,:] = I_.imag

# Or you can do it manually - if you know them; for example:
# (first index contains real/imaginary part)
#I_ = np.zeros([2, model.num_resonances(), model.num_resonances()])
#I_[0, 1, 0] = 2. * 2. / 3.
#I_[0, 0, 0] = 2. * 2.


### DUMP DATA ###
# Check whether background is present
if 'I_background_' in locals():
    data = dict(D = D_, y_data = y_data_, A_cv_data = A_cv_data_, I = I_out_,
                A_v_background_abs2_data = A_v_background_abs2_data_,
                I_background = I_background_)
else:
    data = dict(D = D_, y_data = y_data_, A_cv_data = A_cv_data_, I = I_out_)

stan_rdump(data, MODEL_FOLDER + '/' + args.f_out.name)
print("data_analysis__root_to_dataR.py: Done. Data dumped in {0}.".format(args.f_out.name))


# Make an additional python data file ### DEPRECATED. WILL BE ABANDONED IN THE FUTURE
# (It is useful, if you want to perform analytical, not
# statistical PWA-analysis of the model)

def array_to_string(a):
    """
    Small useful function to convert 2D-arrays to python-readable strings.
    """
    np.set_printoptions(threshold=np.inf)
    res = ' '.join(str(a).split()).replace(' ', ',').replace('\n' , '').replace('[,', '[').replace(',]', ']')
    return 'np.asarray(' + res + ')'

f_py = open(MODEL_FOLDER + '/' + 'amplitude_fitting.data.py', 'w')
f_py.truncate()
#f_py.write('R = ' + str(R) + '\n')
f_py.write('D = ' + str(D_) + '\n')
f_py.write('y_data = ' + array_to_string(y_data_[:,:]) + '\n')
f_py.close()
