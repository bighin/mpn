# Diagrammatic Monte Carlo for electronic correlation in molecules

This repository contains the code developed for the paper G. Bighin, Q.P. Ho, M. Lemeshko, and T. V. Tscherbul, *"Diagrammatic Monte Carlo for electronic correlation in molecules: high-order many-body perturbation theory with low scaling"*. We introduce a novel approach to the calculation of electronic correlation energy in molecules via Diagrammatic Monte Carlo. Please refer to [our preprint](https://arxiv.org/abs/2203.12666) for more details about the physics of the system.

# Requirements

A modern C/C++ compiler (GCC or Clang) and the following libraries:

- [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)
- NCurses
- [psi4](https://psicode.org) is not needed at runtime, but is needed to precompute the electron repulsion integrals.

The following libraries/codes are included in the repository:

- [inih](https://github.com/benhoyt/inih) for parsing .ini files.
- [libprogressbar](https://github.com/doches/progressbar) to display a nice progress bar.

For compiling the code, the CMake build system is needed.

# Usage

The code uses the CMake system, therefore it will be compiled with the commands `mkdir build`, `cd build`, `cmake ..`.

Then go back to the main folder, prepare a .ini file with the details of molecule you want to calculate correlation energies for, the `test.ini` contains an example. Finally run the code (`./build/mpn test.ini`).

# Other information

The folder `psi4` contains script to precalculate the electron repulsion integrals for different molecules. The folder `slurm` contains scripts to run the code on a SLURM cluster.