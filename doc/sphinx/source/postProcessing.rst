.. _postprocessing:

Results Postprocessing
======================

JSON Output Files
-----------------

The majority of output files produced by Phoebe will be in the JSON format. JSON files can be thoughtof as dictionaries. If you process results with Python, for example::

  import json

  # open and parse the file
  with open("the_json_file_name.json","r") as f:
      resultsDict = json.load(f)

  # print the keys of the dictionary to see what the file contains
  print(resultsDict)

This opens a JSON file and prints the keys to the dictionary. The keys indicate which data sets are available in a particular JSON file. We store transport/lifetime data, particle type, band energies, dos, as well as output units in these files.

Plotting Scripts
-----------------

Additionally, in the subfolder ``phoebe/scripts``, we prepared a few simple python scripts that can be used to plot the results in JSON files.

The scripts are self-documented, e.g. you can run::

  ./bands.py -h

to see a description of the script.
In general, these scripts have at least one input argument that must be passed on the command line, including the path to the suitable JSON file.
The scripts generate an output ``pdf`` file at the same location as the input JSON file.

Feel free to use these scripts as a starting point for your own customized plots. If you write an additional script which you think is valuable to other users, consider letting us know and we may add it to the repository.


bands.py
^^^^^^^^^^^^^^^^^^^^^^^^^

Plots the band structure. Usage example::

  ./bands.py path/to/electron_bands.json


transport_coefficients.py
^^^^^^^^^^^^^^^^^^^^^^^^^

Plots the transport coefficients (thermal conductivity, electrical conductivity, etc...) as a function of temperature.
Works for transport coefficients computed with different BTE solvers (EPA, Wannier interpolation, RTA, iterative, ...).
Usage example::

  ./transport_coefficients.py path/to/rta_phonon_thermal_cond.json xx

or::

  ./transport_coefficients.py path/to/epa_onsager_coefficients.json xx


dos.py
^^^^^^^^^^^^^^^^^^^^^^^^^

Plots the density of states. Usage example::

  ./dos.py path/to/electron_dos.json


tau_path.py
^^^^^^^^^^^^^^^^^^^^^^^^^

Overlays the band structure with the particle linewidths, using the output of the lifetimes apps. Usage example::

  ./tau_path.py path_el_relaxation_times.json path_el_bandstructure.json


tau.py
^^^^^^^^^^^^^^^^^^^^^^^^^


Plots lifetimes vs energy. Usage example::

  ./tau.py rta_relaxation_times.json


tau_bandResolved.py
^^^^^^^^^^^^^^^^^^^^^^^^^

Plots lifetimes/linewidths vs energy. It's similar to the script tau.py, but it has a color-code for lifetimes and linewidths from each different band and it only works if no window has been used in the transport code. Note that the n-th band index is defined as the n-th lowest-energy band. Usage example::

  ./tau_bandResolved.py rta_relaxation_times.json



