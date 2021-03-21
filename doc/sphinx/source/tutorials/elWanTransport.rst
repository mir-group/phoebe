Electron Wannier Transport Tutorial
===================================

Synopsis
--------

In this tutorial, we will compute the electrical conductivity and other electronic transport properties of Silicon.
We will use Quantum ESPRESSO to compute the ab-initio electron-phonon coupling on a coarse grid.
Additionally, we will use Wannier90 to compute the maximally localized Wannier functions, which we need to interpolate the electronic band structure and the electron-phonon coupling as well.

The algorithms are described in the Theory section of this manual, and we assume a basic familiarity with Quantum ESPRESSO (see tutorials on phonons on Quantum ESPRESSO's website https://www.quantum-espresso.org/resources/tutorials , and a working knowledge of Wannier90 (although you should still be able to follow this tutorial).


Step 1: Patch Quantum ESPRESSO
------------------------------

As we discussed in the @ref thCOMPILEPHQE page, we need to use a custom modification of Quantum ESPRESSO (which we modified to impose the symmetric properties of the wavefunction).
This will require to compile Quantum ESPRESSO. We provide two options.

* Patch a copy of QE.
  Suppose you have a copy of QE installed at `/oldpath/qe-6.5`.
  Then let's make a copy of that source code (so we don't overwrite your old QE executables) and patch it::

    cd /oldpath
    cp -r qe-6.6 patched-quantum-espresso
    cd patched-quantum-espresso
    make clean
    patch -p1 < /path/to/phoebe/utilities/patch-qe-6.6.txt
    ./configure MPIF90=mpif90 --with-scalapack=yes
    make pw pp ph w90

  where you should modify the file paths to match your installations and you should use the suitable installation parameters for QE that work on your computer.
  If compilation fails, you should consult the QE installation page https://www.quantum-espresso.org/Doc/user_guide/node7.html.

  .. note::
     If you have a different version of QE, you may try to use that patch anyway. We tested it to work for QE-6.5 and QE-6.6. Let us know if you run into problems.

* Alternatively, download our patched version and install it. From an installation folder of your choice, type::

    git clone https://github.com/mir-group/phoebe-quantum-espresso.git
    cd phoebe-quantum-espresso
    git checkout phoebe-qe-6.6
    ./configure MPIF90=mpif90 --with-scalapack=yes
    make pw pp ph w90

where again you should probably customize the installation parameters.





Step 2: pw
----------

First, we need to compute the total energy of the silicon crystal unit cell.
This calculation will create the ground state charge density and wavefunctions that are needed for later.

To run, go to the folder `./example/Silicon-el/qespresso` in the phoebe repository.
The file `scf.in` is the input file for the `pw.x` executable.
The input for a total energy DFT calculation of Quantum ESPRESSO for a silicon crystal, is::

 &control
   calculation = 'scf'
   prefix = 'silicon'
   pseudo_dir = '../../pseudoPotentials/'
   outdir='./out'
   wf_collect = .true.
 /
 &system
   ibrav = 2
   celldm(1) = 10.2
   nat = 2
   ntyp = 1
   ecutwfc = 30.
   nbnd = 12
 /
 &electrons
   conv_thr =  1.0d-14
 /
 ATOMIC_SPECIES
   Si  28.086  Si.pz-vbc.UPF
 ATOMIC_POSITIONS alat
   Si 0.00 0.00 0.00
   Si 0.25 0.25 0.25
 K_POINTS automatic
   4 4 4 0 0 0

A detailed description of all this parameters can be found on Quantum ESPRESSO's website https://www.quantum-espresso.org/Doc/INPUT_PW.html .
The most important parameters to be tweaked and modified in a research project are

* `nbnd`: the number of Kohn-Sham (bands) states to be computed.

  .. note::
     It's important that this parameter is consistent with Wannier90 # of bands.

* `K_POINTS`: parameter controlling the integration mesh of wavevectors in the Brillouin zone. Phonon properties should be converged against this mesh (more wavevectors is better). Tip: `ph.x` calculations are faster when the k-mesh is gamma-centered.

* `ecutwfc`: parameter controlling the number of G-vectors used in the plane-wave expansion of the wavefunction. Phonon frequencies should be checked against this value.

* `conv_thr` this parameter controls total energy convergence. Note that a poorly converged conv_thr may result in poorly converged phonon properties.

* `prefix`: prefix of some output files. Make sure to use a consistent value of prefix throughout your calculations.

* `outdir`: name of the scratch folder. Must be used consistently throughout the run, so to point to the correct files.

This list is obviously not complete, and for your research project you may need to use more functionalities of QE's `pw.x`.

Simply run it as::

    /path/to/patched-quantum-espresso/bin/pw.x -in scf.in > scf.out

after substituting the suitable path to the `pw.x` executable.

.. note::
   We only support the keyword `K_POINTS automatic`.

.. note::
   Use all the bands you need for the Wannierization






Step 3: Phonons and electron-phonon couping
-------------------------------------------

The input file `ph.in` is as follows::

 phonons of Si
 &inputph
  tr2_ph = 1.0d-14
  prefix = "silicon"
  ldisp = .true.
  nq1 = 6, nq2 = 6, nq3 = 6
  outdir = "./out"
  fildyn = "silicon.dyn"
  fildvscf = "silicon.dvscf"
  electron_phonon = "epa"
 /

The values of `nqX` select the Monkhorst-Pack grid of q-points centered at Gamma, for which we will compute the phonon properties.
Here it's important that `prefix` and `outdir` are the same as those used in the `pw.x` calculation of before.
Use a good value of `tr2_ph` (smaller is better, but harder to converge), which (indirectly) checks the convergence of phonon frequencies.

In the input file, we set the flag `electron_phonon = "epa"`.
This will trigger the calculation of the electron-phonon coupling for Phoebe.

Run the code as::

  /path/to/patched-quantum-espresso/bin/ph.x -in ph.in > ph.out

Or in parallel, e.g.::

  mpirun -np 4 /path/to/patched-quantum-espresso/bin/ph.x -npool 4 -in ph.in > ph.out


If the code executes correctly and completely, you should see a number of files called `{fildyn}*`, as many files as the number of irreducible q-points (16 in this case).
On top of that, you should also see several files named as `{prefix}.phoebe.****.dat`, as many as the number of irreducible points.
These files contain the values of the electron-phonon coupling that will be used by Phoebe.

*Current limitations:*

* There are restrictions to the choice of k and q points.
  The `K_POINTS` in `pw.x` must be `automatic`. The `K_POINTS` must be gamma centered.
  And the q-point mesh must be the same as the k-point mesh.

* In the current release, we don't support spin-polarized calculations or spin-orbit calculations. Support for this will come in a later release (we need to implement spin-related symmetries).






Step 4: q2r
-----------

The code ph.x has created the `silicon.dyn*` files, which contain the dynamical matrix at every irreducible q-point.
Now, we run `q2r.x` to Fourier transform the dynamical matrices in the reciprocal space representation to the real space representation, where they represent the interatomic force constants.
The input file `q2r.in` is minimal::

 &input
   fildyn='silicon.dyn',
   flfrc='silicon.fc'
 /

where the first variable must match the path to the dynamical matrices set earlier in `ph.x`, and `flfrc` is the output file with the force constants.

In the working folder `./example/Silicon/qespresso` run the command::

    ./path/to/qe/bin/q2r.x -in q2r.in > q2r.out

If the code run successfully, you should see a new file `silicon.fc`.





Step 5: nscf
------------

We are now moving over to the wannierization of the electronic band structure.
To this aim, we first need to compute the electronic band structure on a full grid of k-points.
You can check that the `nscf.in` file is essentially identical to the `scf.in` file, except that we

* Modified the parameter `calculation = "bands"`, meaning that we will use the charge density computed at step 2 to recompute the wavefunctions.

* Instead of using the keyword `K_POINTS automatic, 6 6 6 0 0 0`, we explicitly write the coordinates of all :math:`6^3` k-points.

To run it, type::

  mpirun -np 4 /path/to/phoebe-quantum-espresso/bin/pw.x -in nscf.in > nscf.out






Step 6: Wannierization
----------------------

Now we can Wannierize the band structure in three steps

First, we run wannier90 in preprocessing mode::

  mpirun -np 4 /path/to/phoebe-quantum-espresso/bin/wannier90.x -pp si


Then we convert data from QE to Wannier90::

  mpirun -np 4 /path/to/phoebe-quantum-espresso/bin/pw2wannier90.x -in pw2wannier90.in > pw2wannier90.out


Finally, the actual wannierization::

  mpirun -np 4 /path/to/phoebe-quantum-espresso/bin/wannier90.x si


The input file of pw2wannier90 is pretty minimal::

 &inputpp
   outdir = './out'
   prefix = 'silicon'
   seedname = 'si'
 /

For your future research project, just make sure that `prefix` and `outdir` are consistent with the `pw.x` calculation above, and that `seedname` is consistent to the name of the wannier90 input file `{seedname}.win`.
The input file of Wannier90 is a bit more involved::

  write_tb = true
  write_u_matrices = true

  num_bands         = 12
  num_wann          = 8
  dis_win_max       = 17.d0
  dis_froz_max      = 6.4d0
  dis_num_iter      = 120
  dis_mix_ratio     = 1.d0

  num_iter          = 500
  num_print_cycles  = 50

  begin unit_cell_cart
  bohr
  -5.1000 0.0000 5.1000
  0.0000 5.1000 5.1000
  -5.1000 5.1000 0.0000
  end unit_cell_cart

  begin atoms_frac
  Si   0.00  0.00   0.00
  Si   0.25  0.25   0.25
  End atoms_frac

  begin projections
  Si : sp3
  end projections

  mp_grid = 6 6 6

  begin kpoints
    0.00000000  0.00000000  0.00000000
    ...
    0.83333333  0.83333333  0.83333333
  end kpoints

Part of this input is just a copy and paste of information coming from the file `nscf.in`.

.. note::
   Additionally:

   you must set the variables::

     write_tb = true
     write_u_matrices = true

   With these, you will write to file the Hamiltonian in the Wannier representation and the rotation matrices :math:`U` that are needed by phoebe.

The variable `num_bands` should match the value of `nbnd` set in `scf.in` and `nscf.in`.

The variable `num_wann` is the number of Wannier functions that we will use in the calculation, and roughly corresponds to the number of bands you can interpolate. The fewer Wannier functions, the faster the calculation.

Finally, in this example we already setup the Wannierization choosing the disentanglement parameters and the orbital projections (that are used as a starting guess of the Wannier orbitals).
The meaning of these quantities is described in the  Wannier90 documentation http://www.wannier.org/support.
Notoriously, this is the hard part of the Wannierization procedure, and every different material may present a new challenge.
The Wannier90 tutorials help you getting started on how to choose these parameters for your research project.



Step 7: QE to Phoebe conversion
-------------------------------

Now that we have all the necessary input files, we can get started with Phoebe.
In this section, we read all the information scattered throughout the files created above and use them to prepare the electron-phonon coupling for the transport calculation.
In detail, we will perform the transformation from the Bloch to the Wannier representation of the electron-phonon coupling.

To do this, let's have a look at the input file `qeToPhoebeWannier.in`::

  appName = "elPhQeToPhoebe"
  elPhInterpolation = "wannier"
  phD2FileName = "silicon.fc"
  electronH0Name = "si_tb.dat"
  wannier90Prefix = "si"
  quantumEspressoPrefix = "silicon"

There are a few parameters to comment:

* :ref:`appName` = `"elPhQeToPhoebe"`:
  here we select the app to postprocess the electron-phonon coupling created by QE.

* :ref:`elPhInterpolation` = `"wannier"`:
  here, we select the postprocessing method that transforms the electron-phonon coupling to the Wannier representation.

* :ref:`phD2FileName` = `"silicon.fc"`

* :ref:`electronH0Name` = `"si_tb.dat"` : this parameter, in the form of `{wannier90seedname}_tb.day` should locate the file created by Wannier90 thanks to the flag `write_tb`. Additionally, there should be present a file called `si_tb_dis.dat` if Wannier90 has disentangled bands.

* :ref:`wannier90Prefix` = `"si"` : should match the `seedname` value of Wannier90, and it is used to locate various `./si.*` files.

* :ref:`quantumEspressoPrefix` = `"silicon"` : this parameter is used to locate and read the files `./silicon.phoebe.*.dat` that have been created by `ph.x`.


There's no other tuning to do, besides making sure that phoebe can locate all these files from this input variables.
To execute the code::

  export OMP_NUM_THREADS=4
  /path/to/phoebe/build/phoebe -in qeToPhoebeWannier.in -out qeToPhoebeWannier.out

and wait until completion.

Note that this calculation can be memory intensive.
For this reason, we recommend to limit/avoid use of MPI parallelization and use a large number of OMP threads (if you compiled the code with OpenMP). (OpenMP facilitates to have multiple threads working on the same memory locations)
MPI parallelization is nevertheless supported also in this code.

After the code completes, you should see an output file called `silicon.phoebe.elph.dat`




Step 8: Electronic Transport with Wannier interpolation
-------------------------------------------------------

Finally, you reached the last step!
Let's see the input file for computing electronic transport properties::

  appName = "electronWannierTransport"
  phD2FileName = "silicon.fc"
  sumRuleD2 = "crystal"
  electronH0Name = "si_tb.dat",
  elphFileName = "silicon.phoebe.elph.dat"

  kMesh = [15,15,15]
  temperatures = [300.]
  dopings = [1.e21]

  smearingMethod = "gaussian"
  smearingWidth = 0.5 eV
  windowType = "population"

  scatteringMatrixInMemory=true
  solverBTE = ["iterative","variational","relaxons"]

There is a number of parameters here:

* :ref:`appName` = `"electronWannierTransport"` : selects the app for computing electronic transport properties with Wannier interpolation.

* :ref:`phD2FileName`: must point to the file `flfrc` created by `q2r.x`, containing the interatomic force constants.

* :ref:`sumRUleD2`: impose translational invariance on the force constants, so that acoustic phonon frequencies go to zero at the Gamma point.

* :ref:`electronH0Name` points to the `si_tb.dat` file created by Wannier90, which contains the electron Hamiltonian in the Wannier representation.

* :ref:`elphFileName` is the path to the file created at step 7 by elPhQeToPhoebe.

* :ref:`kMesh` is the grid used to integrate the Brillouin zone.

  .. note::
     Results must be converged with respect to the `kMesh`

* :ref:`temperatures` in Kelvin, at which we will compute results

* :ref:`dopings` in :math:`cm^{-3}` at which we will compute results. This is only meaningful for semiconductors.

* :ref:`smearingMethod` and :ref:`smearingWidth` sets the algorithm to approximate the Dirac-delta conserving energy. Here we are using the "gaussian" scheme, and the parameter @ref smearingWidth should be converged together with the :ref:`kMesh`. Alternatively, one could use the "adaptiveSmearing" method, which chooses an adaptive width automatically.

* :ref:`windowType` reduces the number of electronic states to those close to the chemical potential. More precisely, selects the electronic states such that :math:`\frac{\partial n}{\partial T} < \delta` and :math:`\frac{\partial n}{\partial \epsilon} < \delta`, where :math:`\delta` is set by :ref:`windowPopulationLimit`. This makes the calculation much faster, as one typically needs just few states close to the chemical potential.

* :ref:`scatteringMatrixInMemory` sets the scattering matrix to be kept in memory, speeding up calculations but costing more memory usage.

* :ref:`solverBTE` selects a number of solvers for the linearized BTE (i.e. solutions beyond the relaxation time approximation, which is always computed)

To run the code, we can simply do::

  export OMP_num_THREADS=4
  /path/to/phoebe/build/phoebe -in electronWannierTransport.in -out ewt.out


Output
------

There are two kinds of output: the output file (in the line above, it's `ewt.out`) and the JSON files.

If you process results with Python, JSON files can be thought as dictionary that are easily read, e.g.::

  import json
  with open("the_json_file_name.json","r") as f:
      resultsDict = json.load(f)

There are several JSON saving all the output, such as the electronic band structure, the electronic lifetimes/linewidths, and the transport properties.

The main output file shows results as well as providing a report of the calculation progress.
The calculation progresses in this way:

* We start by parsing all input files.

* the electronic band structure is computed, and filtered with the window, if needed. In this step, we also compute the Fermi level, chemical potentials and doping concentrations.

* compute the scattering matrix, which is the most time-consuming step.

* Solve the BTE at the relaxation time approximation level, and compute the set of electronic transport coefficients (electrical conductivity, mobility, electronic thermal conductivity and Seebeck coefficient).

<li> Solve the Wigner transport equation at the relaxation time approximation level, and output the transport coefficients.

* Compute the electronic viscosity at the relaxation time approximation level.

* Finally, we start the exact solvers of the linearized BTE. After some time and multiple iterations of the scattering matrix, we compute the transport coefficients. For the "relaxons" solver, we also compute the electronic viscosity obtained by solving the linearized BTE.








Comments
--------

The sections on parallelization discussed for the phonon transport app apply to the electronic transport app as well.

.. note::
   TLDR: :ref:`scatteringMatrixInMemory` = true speeds up calculations but requires a lot of memory (if the code fails the memory allocation, you need to request more HPC resources).
   Using symmetries mitigates the memory usage

   To parallelize your run, set the number of MPI processes equal to the number of nodes, and set the number of OMP threads equal to the number of cores in the node. If applicable, the number of GPUs should match the number of MPI processes.



Here, for simplicity, we are not discussing the convergence tests that need to be done. For a research project on a new untested material, you should make sure to:

* Make sure that phonon frequencies are converged with respect to k-point sampling, q-point sampling and wavefunction cutoff.

* Test the convergence of the Wannier90 bandstructure with respect to the k-point sampling.

* Test the convergence of the electronic transport with respect to ab-initio results, in particular with respect to the k/q-point sampling.

* Check the convergence of the electronic transport results with respect to the parameters @ref kMesh and, if applicable, the @ref smearingWidth.



