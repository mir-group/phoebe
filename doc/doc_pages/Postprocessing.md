@page postproc Results postprocessing
@section postProcScripts Results postprocessing

In the subfolder `phoebe/scripts`, we prepared a few simple python scripts that can be used to plot the results in JSON files.

<ul>
<li> bands.py

Plots the band structure. Usage example:
~~~~~~~~~{.c}
./bands.py path/to/electron_bands.json
~~~~~~~~~

<li> conductivity.py

Plots the thermal conductivity as a function of temperature.
Works for the conductivity computed with different BTE solvers.
Usage example:
~~~~~~~~~{.c}
./conductivity.py path/to/rta_phonon_thermal_cond.json
~~~~~~~~~

<li> dos.py

Plots the density of states. Usage example:
~~~~~~~~~{.c}
./dos.py path/to/electron_dos.json
~~~~~~~~~

<li> tau_path.py

Overlays the band structure with the particle linewidths. Usage example:
~~~~~~~~~{.c}
./tau_path.py path_el_relaxation_times.json path_el_bandstructure.json
~~~~~~~~~

<li> tau.py

Plots lifetimes vs energy. Usage example:
~~~~~~~~~{.c}
./tau.py rta_relaxation_times.json
~~~~~~~~~

</ul>

The scripts are self-documented, e.g. you can run 
~~~~~~~~~~~~~~{.c}
./bands.py -h
~~~~~~~~~~~~~~
to see a description of the script.
In general, these scripts have at least one input argument that must be passed on the command line, and is the path to the suitable JSON file.

The scripts generate an output `pdf` file at the same location of the JSON file.

Feel free to use these scripts as a starting point for your own customized plots and please let us know if you'd like us to include your postprocessing script in Phoebe, if you think it is valuable to other researchers.