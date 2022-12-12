Introduction
============

Phoebe is a software for the ab-initio prediction of transport properties, such as electron-phonon limited conductivity and phonon thermal conductivity.

Currently, we support integration with the ab-initio software suite Quantum ESPRESSO for electron-phonon properties and anharmonic force constants through ShengBTE. Phoebe also accepts anharmonic force constants from Phono3py, which interfaces with most mainstream DFT packages. Support for more coming soon.

Additionally, Phoebe is written in C++ and designed for use on modern HPC infrastructures through hybrid MPI/OpenMP parallelization, distributed memory computation via ScaLAPACK, and support for GPU acceleration using Kokkos.

For more information, or to cite Phoebe, please refer to:

Phoebe: a high-performance framework for solving phonon and electron Boltzmann transport equations.
J. Phys. Mater. 5 035003. (2022). `DOI:10.1088/2515-7639/ac86f6 <https://doi.org/10.1088/2515-7639/ac86f6>`_, `arXiv:2111.14999 <https://arxiv.org/abs/2111.14999>`_.
A. Cepellotti, J. Coulter, A. Johansson, N. S. Fedorova, B. Kozinsky.

.. raw:: html

  <h3>Features:</h3>

* Electron-phonon limited electronic transport coefficients, including the electrical conductivity, Seebeck coefficient, and electronic thermal conductivity

* Lattice thermal conductivity, including effects from 3-phonon scattering, phonon-isotope and phonon-boundary scattering

* Electron-phonon scattering through Wannier interpolation, or EPA approximation

* Calculation of electron and phonon lifetimes (linewidths)

* Electron and phonon viscosities

* Plots of electronic band structure and phonon dispersion relations

* Plots of electron or phonon density of states

This documentation contains a brief description of the formalism used to compute these quantities. We expect that the user is already familiar with ab-initio codes and with the fundamentals of solid state physics. A good introduction to the topic is the book "Electrons and Phonons: The Theory of Transport Phenomena in Solids" by Ziman. In the theory section, we add references to contemporary literature for the latest method/algorithm developments.

This guide will help you through code's installation, calculation of properties, descriptions of the code input files, postprocessing of results, and a theory handbook. A tutorials section will explain how to obtain the transport properties of Silicon. The Developer's documentation section contains a brief overview of the code structure, and is complemented with a detailed description of purpose and methods for each class.
