

### Anharmonic force constants with phono3py
<blockquote>
This is only possible if Phoebe is built using hdf5.
</blockquote>

Alternatively, we also support anharmonic force constants generated from the [phono3py](https://atztogo.github.io/phono3py/) package. This enables phonon transport calculations with Phoebe to be performed using any DFT package with a phono3py interface (VASP, QE, Abinit, etc.)

To calculate the anharmonic force constants from phono3py, you first need to follow the instructions to set up phono3py [here.](https://atztogo.github.io/phono3py/install.html#installation-from-source-code) In many cases, you should be able to do this from a simple conda environment using the following on the command line: 
~~~~~~{.c}
# create a conda environment named phono3py
conda create --name phono3py
# activate the enviroment (use conda activate for newer conda versions)
source activate phono3py
# install phono3py in this environment
conda install -c conda-forge phono3py
~~~~~~
If for any reason this doesn't work, refer to the phono3py documentation linked above for more detailed installation instructions.

The procedure to calculate the force constants with phono3py differs slightly depending on which DFT package you want to use it with. The procedures for interfacing with different calculators are available [here](https://phonopy.github.io/phono3py/interfaces.html#). However, the procedure to follow is similar in all cases. Here, `<input-file-name>` should be a `POSCAR` in the VASP case or a `pw.x` style `pw.in` file in the QE case. `<DFT-package-name>` will be --vasp ( the flag is optional for VASP) and `qe` for Quantum Espresso. For other codes, look at the examples in the above link for input file and package flag information.

Before beginning the calculation, we strongly recommmend you transform whatever unit cell you want to use in the calculation to the primitive cell which will be used internally by phono3py. This avoids the use of flags which signal the need to transform to this primitive cell (see flag [--pa](https://phonopy.github.io/phono3py/command-options.html#pa-primitive-axes-primitive-axes)). If your chosen cell wasn't primitive, this can also reduce the cost of an already expensive calculation.
~~~~~~{.c}
phonopy --<DFT-package-name> --symmetry -c <input-file-name>
~~~~~~
The primitive cell will be written in VASP format in a file named `PPOSCAR`. Put this structure into whatever `<input-file-name>` format you chose. 

Now, we generate the displaced supercells (of the specified dimensions, here 2x2x2) using the command:
~~~~~~{.c}
phono3py --<DFT-package-name> -d --dim="2 2 2" -c <input-file-name>
~~~~~~
and after all the corresponding input files are generated, run a DFT calculation for each of them. 

Once all calculations are finished, collect the force constants from them using a line like:
~~~~~~{.c}
phono3py --<DFT-package-name> --cf3 disp-{00001..nCalculations}/<output-file-name>
~~~~~~
This creates a file named FORCES_FC3, which contains the force constants. To use this information as an input to Phoebe, run the following line to compress this information into two DFT-package independent hdf5 files, `fc2.hdf5` and `fc3.hdf5`, which contain the second and third order force constants, respectively.
~~~~~~{.c}
phono3py --<DFT-package-name> --dim="2 2 2" -c <input-file-name> --sym-fc
~~~~~~

Before proceeding, you should check the quality of the calculation. First, make sure the harmonic phonon bands look appropriate using phono3py. In the directory with the force constants file, make a file named `band.conf` which should contain at a minimum the high symmetry band path in crystal coordinates (with other optional settings [here](https://phonopy.github.io/phonopy/setting-tags.html#band-structure-related-tags)). For silicon, a simple example would be:
~~~~~~{.c}
# save as band.conf
ATOM_NAME = Si
DIM = 3 3 3
BAND = 0.0 0.0 0.0   0.0 0.5 0.5    0.25 0.75 0.5    0.5 0.5 0.5  0.0 0.0 0.0  0.375 0.750 0.375
BAND_LABELS = $\Gamma$ X W L $\Gamma$ K
~~~~~~
Then, run the following lines and check the output plot, named `band.pdf`.
~~~~~~{.c}
phono3py --cfs             # convert the force files to phonopy format  
phonopy --<DFT-package-name> -p -s band.conf -c <input-file-name>     # make a dispersion plot 
~~~~~~
<blockquote>
You should make sure this disperson is converged with respect to DFT convergence parameters (energy cutoff, kpoint mesh, etc) and also the dimension of the supercell provided to phono3py. It is also recommend you check the convergence of the final calculated transport properties with respect to supercell size. 
</blockquote>

If this dispersion looks good, we are now ready to move on to transport calculations using Phoebe. There are four files output by phono3py which we will need: `fc3.hdf5`, `fc2.hdf5`, `phono3py_disp.yaml`, and `disp_fc3.yaml`. These contain all the information we need to go forward, and can be copied into a new directory to run Phoebe if desired. 

Any of the phonon related apps can be run with these files, including the phononBands, phononDos, and lifetime apps. We describe here the use of the transport app here, but the input for other apps will be similar. Below, we provide an example input for the phononTransport app using phono3py inputs. We note that the only two differences are captured by the phD2FileName and phD3FileName inputs. The first input needs to be a path to the directory containing the four files (this applies to all apps, though for those not requiring anharmonic force constants, the phD3FileName variable is not included).
~~~~~~{.c}
appName = "phononTransport"
# below two lines are changed
phD2FileName = "path/to/phono3py/input/directory",
phD3FileName = "fc3.hdf5"

sumRuleD2 = "crystal"
qMesh = [10,10,10]
temperatures = [300.]
smearingMethod = "adaptiveGaussian"
solverBTE = ["variational"]
~~~~~~
Phoebe will then produce the desired transport output in the same way as for the above tutorial using ShengBTE inputs.








