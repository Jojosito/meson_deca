# meson_deca

The project <b>meson_deca</b> is an extension module to CmdStan. Its purpose is to simplify
modelling partial wave analysis of heavy meson physics. We implement such functions as
Breit-Wigner/Flatte dynamical form factors, Blatt-Weisskopf functions and Zemach tensors
and bundle them together to model-dependent PWA amplitudes that describe 3- and 4-body 
meson decays.
  
We introduce only a handful of functions to Stan language;
the core of the module is build from C++ functions that are not
directly callable from Stan (however, they are templated and 
can be exposed to Stan parser, should the need present itself).

The module also contains some small python scripts and modules
that are aimed at debugging the C++ code and plotting/analyzing
the results of the Stan fitting.
  
This code may be interesting to you if you want to:  

 * See an example how to expose your functions to Stan (look at `makefile` and `lib/c_lib/stan_callable`);
 * Fit a function that looks like "f(y,theta) = |A(y) * theta|^2", where "A(y)" and "theta" are complex vectors (look at `lib/c_lib/model.hpp`);
 * Use complex numbers in Stan (look at `lib/c_lib/complex` or `lib/c_lib/stan_callable`);
 * Use some model-dependent PWA functions (templated C++: `lib/c_lib/fct`, Python: `lib/py_lib/fct.py` and `lib/py_lib/res.py`).

### Structure  
Although the main part of the code is, essentially, a Stan-friendly-templated C++ library,
it is not build as a (usual) C++ library. All the C++ code is stored in `lib/c_lib`; 
we do not use forward declarations. We also use the following naming convention: if a
class is defined in "struct_somename.hpp", it is instantiated in "somename.hpp".

### Licensing

Reminder: the core Stan C++ code and CmdStan are licensed under the new BSD.

### Dependencies

* libboost-python-dev (`sudo apt-get install libboost-python-dev` or similar should do the trick);
* PyROOT (you should be able to call `import ROOT` from python shell).

The libboost-python library is used to wrap some C++ functions to a python module (intensely staring
at some python plots can be very useful to find various pesky bugs). The PyROOT package is used by 
some utility scripts that convert Stan output to ROOT trees and vice versa.

### Installation

1. Install the latest CmdStan release (currently supported version: 2.6.2).
  1. Download CmdStan from [https://github.com/stan-dev/cmdstan/releases](https://github.com/stan-dev/cmdstan/releases).
  2. Unpack the tarball.
  3. (Optional) From the CmdStan folder, run 'make' on brief tutorial.

2. Install meson_deca.
  1. From the CmdStan folder, git clone the meson_deca  (type something like `git clone https://github.com/atsipenyuk/meson_deca.git`).
  2. Run `make -s install` from the meson_deca folder.
  3. Add meson_deca/lib/py_lib to your PYTHONPATH. For example, add the following
line to your .bashrc:   
`export PYTHONPATH=$PYTHONPATH:<your path to meson_deca>/meson_deca/lib/py_lib`  


### Usage  

*Main idea.* This module is designed to generate and fit the data according to some distribution
`f_model(y, theta)`.

1. First, you need a `.stan` model that says: "Please generate/fit me some data according to `f_model`".
A template of such a file you can find in `lib/stan_lib`. Go ahead and take a look at
`lib/stan_lib/STAN_data_generator.stan`! (For some other examples, browse `models` folder for `.stan` files.)

2. Each time you will build the `.stan` model, it will look for the function `f_model` in the file
`lib/c_lib/model.hpp`. Make sure that this function is suited to your needs. For details, see the
example below, user guide, or browse the `models/*/backup/model.hpp`.

3. If you think about it, the step 2) above is much more complicated and important than the step 1).
Therefore, when you write your own model, I suggest you proceed as follows.
  1. Adjust the function `f_model` in the file `lib/c_lib/model.hpp`.
  2. Call `./initialize_model FOLDER_NAME` from the `meson_deca` folder. This script will copy the templates of 
the `*.stan` files to `models/FOLDER_NAME` and place a backup copy of `model.hpp` in `models/FOLDER_NAME/backup`.
After that, adjust the `*.stan` files and build them into executables (the script `build.sh` may save you some time).
  3. You can use the executables to sample/fit to your hearts desire. There are some python scripts in `utils/` you may use to convert STAN output to .root files and vice versa.

### Example

(You may want to read `docs/user_guide.pdf` first to get acquainted with the
notation.)

#### Model
This example shows how to use meson_deca to model a decay `D->pi+pi-pi+`
via two toy resonances `f_0(1000)` and `f_0(1200)`. Consider the following
model function:  
`f_model(y,theta) = |theta_1 * A_1(y) + theta_2 * A_2(y) +  theta_3 * A_3(y)|`  
`        = |theta_1 * f_0(1000)(y) + theta_2 * f_0(1200)(y) + theta_3 * 1|^2.`  

This model is implemented in the file `models/d_to_3pi_model_dep/backup/model.hpp`. To use the model, copy this
file to `lib/c_lib/model.hpp`:  
` ../meson_deca $ cp models/d_to_3pi_model_dep/backup/model.hpp lib/c_lib`  
The file `lib/c_lib/model.hpp` is the main file you have to adjust when you want to make your own model. There are several key points:

1) Since we are considering a 3-body-decay, the variable `y` is two-dimensional: `y = (m2_ab, m2_bc)`. In `model.hpp`, this is fixed via the line `NUM_VAR=2`.  

2) There are 3 resonances: `A_1`, `A_2`, and `A_3`. In `model.hpp`, the number
of resonances is fixed via the line `NUM_RES=3`.  

(The 3rd resonance is just a dummy: in the end, we shall set theta_3 = 0. 
We keep this resonance in the example code just so `NUM_RES` and `NUM_VAR` are
not equal (because they are, in general, completely independent from each 
other)).  

3) In `model.hpp`, the resonances are fixed in the function `A_c`; each resonance returns a complex number.  

4) In `model.hpp`, the resonances are bundled together in the function `A_cv`, which returns the complex vector `(A_1, A_2, ..., A_NUM_RES)`.  

#### Goal
First, we want to fix `theta = (theta_1_init, theta_2_init, 0)` and generate 
10000 events `{y_i}` according to our model. Then, we want to use the 
generated data to sample the complex parameter `theta_2` (the parameter 
`theta_1 = theta_1_init` remains fixed as reference parameter). The result 
of this sampling should be sharpely peaked around the value `theta_2_init`.

#### Code

Copy pre-build STAN files into a new folder bw2_example and 
wrap model to a python module. (Some warnings may be casted, but
no errors should occur. NOTE: the code is wrapped to a python
module using `clang++`. If you would like to use `g++`, adjust the file
`meson_deca/lib/c_lib/py_wrapper/setup.py`)  
  
`../meson_deca $ ./initialize_model.sh bw2_example`  
`../meson_deca $ cd models/bw2_example`  
  
Build the STAN files into executables.  
`../bw2_example $ ./../../build.sh`  

Pre-calculate the normalization integral I (cf. docs/user_guide.pdf);
basically, calculate `\int A_i(y) A_j(y) dy`. The variable `y` is a 2-dim.
vector, so we need to pass integration boundaries for `y.1` and `y.2`.  
`../bw2_example $ ./../../utils/calculate_normalization_integral.py 0 3 0 3`  
  
Generate 10000 events (STAN puts them into a `generated_data.csv` file);
plot the results; convert .csv file to a .root file with trees y.1, y.2;
convert .root tree to the `STAN_amplitude_fitting.data.R` file, which 
contains `(A_1(y), ... A_3(y))` for 10000 events y.  
`../bw2_example $ ./../../generate.sh 10000`  
  
You can look at the plotted data:  
`../bw2_example $ evince generated_data.pdf`  
  
Use the data file `STAN_amplitude_fitting.data.R` to sample 3 chains, each
with 1000 values of theta_2:  
`../bw2_example $ ./../../fit.sh 1000`  
  
After the simulation is finished, plot the sampled theta values:  
`../bw2_example $ ./../../plot_fit.sh`  
`../bw2_example $ evince output.pdf`  




