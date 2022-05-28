.. _phononTransport:

Phonon Transport Tutorial
=========================

Synopsis
--------

While in many semiconducting and insulating systems phonon transport is dominated by phonon-phonon scattering, in metals and some systems at low temperature, it's necessary to include both the phonon-phonon and phonon-electron scattering in the calculation of lattice thermal conductivity. Additionally, in metallic systems the electronic part of the thermal conductivity often provides the largest contribution to the thermal conductivity. Because Phoebe offers both electron-phonon and phonon-phonon interactions, we can calculate the full thermal conductivity :math:`\kappa` of a metal using phonon-phonon and phonon-electron scattering for :math:`\kappa_{ph}` as well as electron-phonon scattering to calculate :math:`\kappa_{el}`.

This is no small task, as the calculation requires the input files from both the phonon-phonon and electron-phonon calculations. Before starting, we ask that you follow the tutorials for :ref:`phononTransport` (steps 1-3) and :ref:`elPhTransport` (steps 1-7), using the input files `here <phoebe site>`_ for copper.

While we note that these tutorial files will result in unconverged calculations, they will serve the purpose of this example. Both sets of calculations need to be converged as described in their respective tutorials for real use. 

To begin, we will require these input files from each of the coupling calculations: 

**Electron-Phonon**: 
* `cu.fc`
* `cu_tb.dat`
* `cu.phoebe.elph.dat`

**Phonon-Phonon**: 
* `fc2.hdf5`
* `fc3.hdf5`
* `phono3py_disp.yaml` 
(though ShengBTE files could also be used).


Step 1: Electronic Thermal Conductivity Calculation
----------------------------------------------------

First we do the straightforward step -- we need to calculate :math:`\kappa_{el}`. This can be done as in the :ref:`calculateElTransport` step of the electron Wannier transport tutorial. All of the variables we will use are described in that tutorial, so be sure to look over what is written there. 

Then, run the following Phoebe input file:: 

  appName = "electronWannierTransport"
  phFC2FileName = "cu.fc"
  sumRuleFC2 = "crystal"
  electronH0Name = "cu_tb.dat",
  elphFileName = "cu.phoebe.elph.dat"

  kMesh = [20,20,20]
  temperatures = [100.,200.,300.,400.]
  chemicalPotentials = [FERMI ENERGY]

  smearingMethod = "adaptiveGaussian"
  windowType = "population"

  useSymmetries = true
  scatteringMatrixInMemory = false

Here, we do a calculation across a few temperatures, and therefore, for the sake of doing a quick calculation, we use the RTA transport method. This means we don't need to store the full scattering matrix in memory and we can run them all together to save time. In a production quality calculation, it would also be important to be sure each temperature is converged for the selected k-point mesh, as lower temperatuers will require more k-points to converge. 

This should be run just as in the other tutorial::

  export OMP_NUM_THREADS=4
  /path/to/phoebe/build/phoebe -in electronWannierTransport.in > electronTransport.out

Step 2: Lattice Thermal Conducivity Calculation 
------------------------------------------------

Now, we have to calculate the lattice thermal conductivity using both phonon-electron and phonon-phonon scattering. This corresponds to :ref:`calculatePhononTransport`, step 4 of the phonon transport tutorial, but now also uses electron-phonon information. 

For this, we run a Phoebe calculation using the below input file::

  appName = "phononTransport"

  # specify the paths to phonon-phonon input
  phFC2FileName = "fc2.hdf5"
  phFC3FileName = "fc3.hdf5"
  phonopyDispFileName = "phono3py_disp.yaml"

  # tell phoebe to read in the electron-phonon information
  usePhElScattering = true

  # path to electron-phonon coupling inputs 
  phFC2FileName = "cu.fc"  # TODO need to add an extra item for this....? 
  electronH0Name = "cu_tb.dat",
  elphFileName = "cu.phoebe.elph.dat"

  temperatures = [100.,200.,300.,400.]
  chemicalPotentials = [FERMI ENERGY]
  sumRuleFC2 = "crystal"
  qMesh = [15,15,15]
  smearingMethod = "adaptiveGaussian"

  useSymmetries = true
  scatteringMatrixInMemory = false

Again, this should be run just as in the other tutorial::

  export OMP_NUM_THREADS=4
  /path/to/phoebe/build/phoebe -in phononTransport.in > phononTransport.out


Step 3: Post-Process the Outputs
------------------------------------------------

From these two calculations, we'll need the `solver_onsager_coefficients.json` from the electronic calculation, and the `rta_phonon_thermal_cond.json` files. 
Below, we plot the output of these calculations together using the following simple python script:: 

  # TODO python script


# TODO stick in the plot 





