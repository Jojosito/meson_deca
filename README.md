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
  1a. Download CmdStan from [https://github.com/stan-dev/cmdstan/releases](https://github.com/stan-dev/cmdstan/releases).
  1b. Unpack the tarball.
  1c. From the CmdStan folder, run 'make' on brief tutorial.

2a. Clone the meson_deca into the CmdStan folder (type something like `git clone https://github.com/atsipenyuk/meson_deca.git`).

2b. Run `make install` from the meson_deca folder.
