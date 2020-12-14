Introduction
============

Phoebe is a software for the characterization of transport properties from ab-initio, such as electron-phonon limited conductivity and phonon thermal conductivity. This C++ code supports GPU acceleration for modern HPC infrastructures.

Currently, we fully support integration with the ab-initio software suite Quantum ESPRESSO, with support for more coming soon. Features include:

* Electron transport coefficients (electron conductivity, Seebeck coefficient, electronic thermal conductivity, limited by electron-phonon scattering;

* Lattice thermal conductivity, including effects from 3-phonon scattering, phonon-isotope and phonon-boundary scattering;

* Calculation of electron or phonon lifetimes (linewidths);

* Electron and phonon viscosity;

* Plot of electronic band structure and phonon dispersion relations;

* Plot of electron or phonon density of states;

* Electron-phonon scattering through Wannier interpolation, or EPA approximation.

This documentation contains a brief description of the formalism used to compute these quantities. We expect however that the user is already familiar with ab-initio codes and with the fundamentals of solid state physics. A good introduction to the topic is the book "Electrons and Phonons: The Theory of Transport Phenomena in Solids" by Ziman. In the theory section, we add references to contemporary literature for the latest method/algorithms developments.

This guide will help you through code's installation, calculation of properties, descriptions of the code input files, postprocessing of results, and a theory handbook. A tutorials section will explain how to obtain the transport properties of Silicon. The Developer's documentation section contains a brief overview of the code structure, and is complemented with a detailed description of purpose and methods for each class.
