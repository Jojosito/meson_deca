# meson_deca

<b>meson_deca</b> is an extension module to Stan. Its purpose is to simplify
modelling partial wave analysis of heavy meson physics.

### Licensing

Reminder: the core Stan C++ code and CmdStan are licensed under new BSD.

### Dependencies

* libboost-python-dev
* PyROOT

### Installation

1. Install the latest CmdStan release (current support: 2.6.2).
  1. Download CmdStan from [https://github.com/stan-dev/cmdstan/releases](https://github.com/stan-dev/cmdstan/releases).
  2. Unpack the tarball.
  3. (Optional) From the CmdStan folder, run 'make' on brief tutorial.

2. Install meson_deca.
  1. From the CmdStan folder, git clone the meson_deca  (type something like `git clone https://github.com/atsipenyuk/meson_deca.git`).
  2. Run `make -s install` from the meson_deca folder.

### Example

(You may want to read `docs/user_guide.pdf` first to get acquainted with the
notation.)

#### Model
This example shows how to use meson_deca to model a decay D->pi+pi-pi+
via two toy resonances f_0(1000) and f_0(1200). Consider the following
model function:  
`f_model(y,theta) = |theta_1 * A_1(y) + theta_2 * A_2(y) +  theta_3 * A_3(y)|`  
`        = |theta_1 * f_0(1000)(y) + theta_2 * f_0(1200)(y) + theta_3 * 1|^2.`  

You can check how this model is implemented by looking at the file `lib/c_lib/model.hpp` - this is the main file you have to adjust when you want to make your own model. There are several key points:

1) Since we are considering a 3-body-decay, the variable y is two-dimensional:
y = (m2_ab, m2_bc). In `model.hpp`, this is fixed via the line `NUM_VAR=2`.  

2) There are 3 resonances: A_1, A_2, and A_3. In `model.hpp`, the number
of resonances is fixed via the line `NUM_RES=3`.  

(The 3rd resonance is just a dummy: in the end, we shall set theta_3 = 0. 
We keep this resonance in the example code just so NUM_RES and NUM_VAR are 
not equal (because they are, in general, completely independent from each 
other).  

3) The resonances are fixed in the function A_c; each resonance returns
a complex number.  

4) The resonances are bundled together in the function A_cv, which returns
the complex vector (A_1, A_2, ..., A_NUM_RES).  

#### Goal
First, we want to fix theta = (theta_1_init, theta_2_init, 0) and generate 
10000 events {y_i} according to our model. Then, we want to use the 
generated data to sample the complex parameter theta_2 (the parameter 
theta_1 = theta_1_init remains fixed as reference parameter). The result 
of this sampling should be sharpely peaked around the value theta_2_init.

#### Code

Copy pre-build STAN files into a new folder bw2_example and 
wrap model to a python module. (Some warnings may be casted; but
no errors should occur.)  
  
../meson_deca $ ./initialize_model bw2_example  
../meson_deca $ cd models/bw2_example  
  
Build the STAN files into executables.  
../bw2_example $ ./../../build.sh  

Pre-calculate the normalization integral I (cf. docs/user_guide.pdf);
basically, calculate `\int A_i(y) A_j(y) dy`. The variable `y` is a 2-dim.
vector, so we need to pass integration boundaries for `y.1` and `y.2`.  
../bw2_example $ ./../../utils/calculate_normalization_integral.py 0 3 0 3  
  
Generate 10000 events (STAN puts them into a 'generated_data.csv' file);
plot the results; convert .csv file to a .root file with trees y.1, y.2;
convert .root tree to the `STAN_amplitude_fitting.data.R` file, which 
contains (A_1(y), ... A_3(y)) for 10000 events y.  
../bw2_example $ ./../../generate.sh 10000  
  
You can look at the plotted data:  
../bw2_example $ evince generated_data.pdf  
  
Use the data file `STAN_amplitude_fitting.data.R` to sample 3 chains, each
with 1000 values of theta_2:  
../bw2_example $ ./../../fit.sh 1000  
  
After the simulation is finished, plot the sampled theta values:  
../bw2_example $ ./../../plot_fit.sh  
../bw2_example $ evince output.pdf  




