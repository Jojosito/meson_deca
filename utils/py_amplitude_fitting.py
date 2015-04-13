#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import os
import sys

MODEL_FOLDER = os.getcwdu()
sys.path.insert(1, MODEL_FOLDER)
import model # Contains A_r - model-dependent PWA resonance function

sys.path.insert(1, "../../lib/py_lib")
import convert # Translates A_r results to usable form


def Atheta(y, theta):
    # y is a float vector with model.num_variables() entries
    # theta is a complex vector with model.num_resonances() entries
    return np.dot(convert.ComplexVectorForm(model.A_r(model.num_variables(),y)), theta)


def H(y, theta):
    return np.square(np.abs(Atheta(y, theta)))


# Theta prior - assuming uniform distribution
# theta_x will be built-in in theta
theta_x_min = -2.
theta_x_max = 2.
theta_x_list = np.linspace(theta_x_min, theta_x_max, 100)
theta_list = [[1., theta, 0.] for theta in theta_x_list]

# x range (one-dimensional x will be built-in in y)
x_min = -1.
x_max = 1.
x_list = np.linspace(x_min, x_max, 80)


execfile('normalization_integral.py')
def Norm(theta):
    #return scipy.integrate.quad(lambda y: f(y,theta), y_min, y_max)[0]
    return np.dot(theta, np.dot(I_,theta))

def f_norm(y, theta):
    return H(y, theta) / Norm(theta)


# Posterior
def Pr(theta, y):
    return f_norm(y, theta)

y_list = [[x,0] for x in x_list]

def Pr_list(theta, y_list):
    return np.sum([np.log(Pr(theta, y_list[:][n])) for n in range(len(y_list))])


def Pr_sigma(y_list):
    return np.asarray([Pr_list(theta, y_list) for theta in theta_list])


plt.clf()
z_0 = np.asarray([f_norm(y, [1.,0.,0.]) for y in y_list]) 
z_m2 = np.asarray([f_norm(y, [1.,-2.,0]) for y in y_list])
z_2 = np.asarray([f_norm(y, [1.,2.,0.]) for y in y_list])
z_m1 = np.asarray([f_norm(y, [1.,-1.,0.]) for y in y_list])
z_1 = np.asarray([f_norm(y, [1.,1.,0.]) for y in y_list])
plt.plot(x_list, z_m2, label = 'theta = -2.')
plt.plot(x_list, z_m1, label = 'theta = -1..')
plt.plot(x_list, z_0, label = 'theta = 0.')
plt.plot(x_list, z_1, label = 'theta = 2.')
plt.plot(x_list, z_2, label = 'theta = 1.')
plt.legend()
plt.savefig('model.pdf')
plt.clf()


if 'y_data' in locals():
    del(y_data)

plt.figure()
execfile(MODEL_FOLDER + '/' + 'amplitude_fitting.data.py')
logp = Pr_sigma(y_data.tolist())
plt.plot(theta_list, logp)
plt.savefig(MODEL_FOLDER + '/' + 'log_likelihood.pdf')

plt.clf()
plt.figure()
plt.plot(theta_list, np.exp(10*logp/abs(max(logp))))
plt.savefig(MODEL_FOLDER + '/' + 'likelihood.pdf')










