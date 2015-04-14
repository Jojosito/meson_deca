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

1c.(Optional) From the CmdStan folder, run 'make' on brief tutorial.

2a. From the CmdStan folder, git clone the meson_deca  (type something like `git clone https://github.com/atsipenyuk/meson_deca.git`).

2b. Run `make -s install` from the meson_deca folder.

### Example

