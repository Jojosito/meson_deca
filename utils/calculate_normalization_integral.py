#!/usr/bin/env python
# calculate_normalization_integral.py

# NAME
#    calculate_normalization_integral.py - make a file containing I[i,j]
#
# SYNOPSIS
#    ./calculate_normalization_integral.py [OPTIONS...]
#
# DESCRIPTION
#    calculate_normalization_integral.py calculates the integrals
#    necessary to perform the correct normalization of the
#    log-likelihood function during STAN_amplitude_fitting
#    sampling.
#
#    This means the following: the function model.A_r (python-wrapped
#    analogon of the STAN function A_cv) contains R =
#    model.num_resonances() amplitudes. This script calculates
#
#        I[i,j] = \int conj(A_r(i)) A_r(j) \dy,
#
#    where y is the model variable - for examle, the invariant mass of
#    the decay. The integral is calculated using Monte-Carlo
#    integration and is stored in 'normalization_integral.py'. The
#    integration bounds may be specified in [OPTIONS] - just pass the
#    arguments in the form y1_min, ... , yR_min, y1max, ... , yR_max.
#
# CAVEAT. This script MUST be called from the folder containing the
#    module model.so corresponding to the described model.

import argparse
import numpy as np
import os
import sys

# Next line is bad (rework relative path)
sys.path.insert(1, os.getcwdu() + '/' + "../../lib/py_lib")
import mcint   # Monte Carlo integration
import convert # Translates A_r results to usable form
import save    # Convert an array to string

# Import PWA data from current model
sys.path.insert(1, os.getcwdu())
import model  # Function to be integrated

# Parse the arguments
parser = argparse.ArgumentParser(description='Script to calculate the integral leading to the Normalization intergral.')

# Parse arguments - two boundaries for each resonance!
# y1min, y1max, ... yRmin, yRmax
N = model.num_variables()
R = model.num_resonances()
parser.add_argument('bounds',
                     nargs=2*N,
                     type=float,
                     help="integration bounds.")

args = parser.parse_args()

I = np.zeros([R,R], dtype=complex)

# Complex vector containing all PWA amplitudes
def func(y):
    return convert.ComplexVectorForm(model.A_cv(N, y))


# Repack the integration bounds
bounds = [[args.bounds[2*n], args.bounds[2*n+1]] for n in range(N)]

# Calculate the integral matrix
I = mcint.integral_A(func, bounds, N=1000000)

f_py = open('normalization_integral.py', 'w')
f_py.write('I_ = ' + save.array_to_string(I[0]) + '\n')
f_py.close()
