#!/usr/bin/env python

import os
from distutils.core import setup
from distutils.extension import Extension

# Switch between g++/clang++ by altering
# the next two lines
os.environ["CC"] = "clang++"
os.environ["CXX"] = "clang++"

# set cmdstan directory (as the first directory containing '/meson_deca/')
cmdstan_path = os.getcwd()
cmdstan_path = cmdstan_path[:cmdstan_path.index('/meson_deca/')]

setup(name="Model_Dep_Functions",
      ext_modules=[
          Extension("model", ["model.cpp"],
          libraries = ["boost_python"],
          include_dirs=[cmdstan_path + "/stan/lib/boost_1.55.0",
                        cmdstan_path + "/stan/lib/eigen_3.2.4",
                        cmdstan_path + "/stan/src",
                        cmdstan_path])
      ])
