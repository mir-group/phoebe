Phonon Transport Tutorial
=========================

Synopsis
--------

In this tutorial, we want to compute the lattice thermal conductivity of Silicon.
We will use Quantum ESPRESSO as the code that provides DFT parameters.
As described in the :ref:`Theory` section of this manual, we need Quantum ESPRESSO to compute the phonon dynamical matrix and the third derivative of the total energy with respect to ionic displacements.

Note that we assume the reader to be familiar with Quantum ESPRESSO.
Several tutorials can be found on Quantum ESPRESSO's website https://www.quantum-espresso.org/resources/tutorials , which cover more complicated DFT calculations than that described here.


Step 1: pw.x, total energy calculation
--------------------------------------

First, we need to compute the total energy of the silicon crystal unit cell.
This calculation will create the ground state charge density and wavefunctions that are needed for later.

To run, go to the folder `./example/Silicon/qespresso` in the phoebe repository.
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

A detailed description of all this parameters can be found on Quantum ESPRESSO's website https://www.quantum-espresso.org/Doc/INPUT_PW.html.
The most important parameters to be tweaked and modified in a research project are

* `K_POINTS`: parameter controlling the integration mesh of wavevectors in the Brillouin zone. Phonon properties should be converged against this mesh (more wavevectors is better). Tip: `ph.x` calculations are faster when the k-mesh is gamma-centered.

* `ecutwfc`: parameter controlling the number of G-vectors used in the plane-wave expansion of the wavefunction. Phonon frequencies should be checked against this value.

* `conv_thr` </li> this parameter controls total energy convergence. Note that a poorly converged conv_thr may result in poorly converged phonon properties.

* `prefix`: prefix of some output files. Make sure to use a consistent value of prefix throughout your calculations.

* `outdir`: name of the scratch folder. Must be used consistently throughout the run, so to point to the correct files.

This list is obviously not complete, and for your research project you may need to use more functionalities of QE's `pw.x`.

Simply run it as::

  /path/to/qe/pw.x -in scf.in > scf.out

after substituting the suitable path to the `pw.x` executable.





Step 2: ph.x, phonon dynamical matrices
---------------------------------------

The input file `ph.in` is as follows::

 phonons of Si
 &inputph
  tr2_ph = 1.0d-14
  prefix = "silicon"
  ldisp = .true.
  nq1 = 4
  nq2 = 4
  nq3 = 4
  outdir = "./out"
  fildyn = "silicon.dyn"
 /

The values of `nq*` select the Monkhorst-Pack grid of q-points centered at Gamma, for which we will compute the phonon properties.
Here it's important that `prefix` and `outdir` are the same as those used in the `pw.x` calculation of before.
Use a good value of `tr2_ph` (smaller is better, but harder to converge), which (indirectly) checks the convergence of phonon frequencies.

Run the code as::

    /path/to/qe/bin/ph.x -in ph.in > ph.out

If the code executes correctly and completely, you should see a number of files called `{fildyn}*`, as many files as the number of irreducible q-points.

.. note::
   We strongly recommend to parallelize this calculation, and to use `npool`. This parameter controls the parallelization over k-vector. It gives an almost linear speedup, but is bound by the number of k-points. For example::

     mpirun -np 4 /path/to/qe/bin/ph.x -npool 4 -in ph.in > ph.out


.. note::
   This is a simple example built such that the parameters used here yield reasonable results.
   In general, we strongly recommend to test convergence of the phonon frequencies with respect to the k-point mesh, the q-point mesh, the wavefunction cutoff and the `tr2_ph` parameter.
   We also recommend to use a small `conv_thr`.






Step 3: q2r.x, harmonic force constants
---------------------------------------

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




Step 4: anharmonic force constants
----------------------------------

In this section, we want to use a finite-displacement approach to computing the matrix of third derivatives of the total energy with respect to ionic displacements.
To this aim, we will be using Quantum ESPRESSO to compute energies/forces, and thirdorder.py to generate a pattern of displacements in a supercell of the crystal.

* Download thirdorder.py from here http://www.shengbte.org/downloads

* Untar the file and cd into the `./thirdorder` directory that has been just created

* Modify the source code in the following way.
  Modify line 559 of file thirdorder_core.c, from `#include "spglib/spglib.h"` to `#include "spglib.h"`.
  In file setup.py, set line 10 as `INCLUDE_DIRS = ["/your/path/to/phoebe/build/spglib_src/src"]` and line 13 as `LIBRARY_DIRS = ["/your/path/to/phoebe/build/spglib_build"]`.

* Open a terminal in the `thirdorder` directory and type::

    ./compile.sh

  If everything works, you should find a `*.so` file in the subdirectories of `./build`.

* Let's go back to the qespresso directory `/path/to/phoebe/example/Silicon/qespresso`.
  Let's check the file `supercell_template.in`.
  The content should look as::

    &control
      calculation = 'scf'
      restart_mode='from_scratch',
      prefix='silicon',
      tstress = .true.
      tprnfor = .true.,
      pseudo_dir = '../../pseudoPotentials/',
      outdir='./out',
    /
    &system
      ibrav = 0
      nat = ##NATOMS##
      ntyp = 1,
      ecutwfc = 30.
    /
    &electrons
      conv_thr =  1.0d-12
    /
    ATOMIC_SPECIES
      Si  28.086  Si.pz-vbc.UPF
    ##COORDINATES##

    ##CELL##
    K_POINTS gamma

  As you can notice, the file is the same as `scf.in`, but we modified a few things:

   * we set `tstress` and `tprnfor` to true.

   * we removed `celldm` (and you should remove `alat`, if used)

   * we set `ibrav=0`

   * we set a tag in place of the number of atoms `nat`.

   * Removed Cell and Coordinates cards and replaced them with tags

   * Modified the k-points, as the k-point density should decrease like the size of the supercell we will set up. In this case, we initially set a k-point mesh of 4x4x4 points, but we will set up a supercell of size 4x4x4 and thus the new supercell k-point mesh is 1x1x1.



.. note::
   If you use the `K_POINTS gamma` keyword, make sure you don't use the patched version of QE modified for the electron-phonon coupling, or use it with `K_POINTS automatic`.


* Now, we generate the displacements on the supercell that are needed to compute the third-order force constants.
  From the phoebe example directory, run in the terminal::

    ln -s /your/path/to/thirdorder_espresso.py .
    python3 thirdorder_espresso.py scf.in sow 4 4 4 -3 supercell_template.in

  In the first command, we link the script provided by `thirdorder`, please modify it to match the correct path.
  Next, you can see the script takes 7 parameters.

     * First, the QE input for the unit cell.

     * Next, `sow` means we generate the supercells

     * 4 4 4 is the three parameters indicating the 4x4x4 supercell size

     * -3 indicates that we only include interactions up to the third nearest neighbor.

     * Finally, we pass the path to the supercell template discussed above

  This script will create a lot of input files, potentially, up to the cube of the number of atoms in the supercell, therefore choose an appropriate number of nearest neighbors (by converging the thermal conductivity)!

* Now, it's time to run all of these supercell calculations!
  For example, you can do this by typing in the terminal::

    for f in DISP.supercell_template.in.*; do
      mpirun -np 4 pw.x -in $f > $f.out
    done

  This step may take a while...

* Finally, we postprocess all these forces by typing::

    find . -name 'DISP.supercell_template.in.*out' | sort -n | python3 thirdorder_espresso.py scf.in reap 4 4 4 -3

  Note here that you should use the same parameters (here, 4 4 4 -3) used for generating the supercell displacements.
  If everything goes well, you should see a new file called `FORCE_CONSTANTS_3RD` with the desired output.

Congratulations! You computed the ab-initio matrix of third order force constants.







Step 5: Phoebe, phonon transport
--------------------------------

The typical input file looks like this::

  appName = "phononTransport"
  phD2FileName = "./qe-phonons/silicon.fc",
  phD3FileName = "./qe-ph-anharmonic/FORCE_CONSTANTS_3RD"
  sumRuleD2 = "crystal"
  qMesh = [10,10,10]
  temperatures = [300.]
  smearingMethod = "adaptiveGaussian"
  solverBTE = ["variational"]


Let's go through this parameters one by one:

* :ref:`appName` = `"phononTransport"` triggers the calculation of phonon transport properties

* :ref:`phD2FileName` must point to the `flfrc` file produced by `q2r.x`

* :ref:`phD3FileName` must point to the file of third derivatives

* :ref:`sumRuleD2` allows us to re-enforce the translational-invariance of force constants, that is broken by numerical errors. After imposing this sum rule, acoustic phonon frequencies to go to zero at the gamma point.

* :ref:`qMesh` is the size of the grid of wavevectors used to integrate the Brillouin zone. Note that the value used here is very unconverged, so that the example can finish in a short amount of time.

  .. note::
     Results must be converged against values of `qMesh`!

* :ref:`temperatures` sets the temperature in Kelvin

* :ref:`smearingMethod` sets the algorithm to approximate the Dirac-delta conserving energy. Using the "adaptiveGaussian" scheme is particular convenient as the width of the gaussian is automatically adjusted. With "gaussian" scheme instead, you should converge the :ref: `smearingWidth` parameter together with the :ref:`qMesh`.

* :ref:`solverBTE` Selects the algorithm to solve the linearized Boltzmann equation. If not specified, we only compute results within the relaxation time approximation. Here, we are using the variational solver to find the solution to the BTE.

With this input, we can compute the phonon contribution to thermal conductivity of silicon.

.. note::
   By default, isotopic scattering at natural abundances is included in the scattering matrix. To disable or modify it, check the parameters :ref:`withIsotopeScattering` and :ref:`massVariance`.


.. note::
   In several studies you may want to include boundary scattering. To include it, use the parameter :ref:`boundaryLength`.





Output
------

Here is what the code is doing:

* parsing input files

* Computing the phonon band structure (energies, eigenvectors and velocities)

* Computes the scattering matrix (this takes place whenever you see a block like this one::

    Started computing scattering matrix with 64 q-points.
    2020-10-30, 09:15:02 |   1% |  1 / 64
    2020-10-30, 09:15:02 |   4% |  3 / 64
    2020-10-30, 09:15:02 |   9% |  6 / 64 | remaining: 6.62e-01 s.
    ......
    2020-10-31, 09:15:03 | 100% | 64 / 64 | remaining: 2.50e-02 s.
    Elapsed time: 0.81 s.

  where, for your convenience, we try to estimate the time to completion.

* Thermal conductivity, BTE theory, estimated within the relaxation time approximation.

* Wigner Thermal conductivity, obtained including off-diagonal contributions of the flux operator, estimated within the relaxation time approximation.

* Thermal viscosity tensor within the relaxation time approximation.

* Lattice contribution to specific heat (at constant volume)

* Optional: if you selected an exact solver, you will see additional output, which includes the thermal conductivity obtained by solving the full linearized BTE (including off-diagonal matrix elements of the scattering operator).

* Optional: if you use the relaxon solver, you will also see the thermal viscosity obtained by solving the BTE exactly.

Note also that the code write results in a variety of JSON files, for ease of use.
If, for example, you use Python for result postprocessing, you can load them as::

  import json
  with open("rta_phonon_thermal_cond.json") as f:
    a=json.load(f)

After this lines, the JSON is loaded in the variable `a` as a dictionary and is ready to be postprocessed.




Tradeoffs between speed and memory
----------------------------------

There's a parameter :ref:`scatteringMatrixInMemory` that you need to consider.
If we set this parameter to true, we store the scattering matrix in memory.
If false, we only compute the action of the scattering matrix, without ever storing all of it in memory.

There is no `best` choice here, rather, you should decide what's best for your case and decide which tradeoff works best for you.

* Option 1: :ref:`scatteringMatrixInMemory` = true. The scattering matrix occupies :math:`16 (3 N_{atoms} N_{q-points})^2 / 1024^3` Gigabytes, if no window is used. This number can be pretty large (even Terabytes), and you should make sure that your HPC allocation has enough memory for storing this large matrix. Given the size, we only allow you to run the code with a single temperature.

  In exchange, iterative or variational solvers of the BTE are extremely cheap, and the cost of your simulation is basically the cost of constructing the scattering matrix. Moreover, this allows you to run :ref:`solverBTE` = "relaxons" type of BTE solver.

* Option 2: :ref:`scatteringMatrixInMemory` = false. The memory footprint is much lighter (the square root of before), so that the same calculation can be run on fewer CPUs. You can compute the thermal conductivity for multiple temperatures in the same run. The calculation of properties within the relaxation time approximation is as expensive as above (if this is what you care about, definitely use less memory).

  In exchange, iterative or variational BTE solvers are much slower. In fact, at each iteration you need to recompute the scattering matrix.
The cost of the calculation therefore grows linearly with the number of iterations of the iterative solver (which may be significant).
You also cannot diagonalize the scattering matrix with :ref:`solverBTE` = "relaxons".





Low temperature thermal conductivity
------------------------------------

At low temperatures, only phonons with small energies are thermally excited and most states are empty.
However, if we use the same input file as the one above, we are sampling all phonon states.
As a result, we end up spending a lot of time computing phonon states that don't contribute to transport.

To address this, we have the parameters :ref:`windowType`, :ref:`windowEnergyLimit`, and :ref:`windowPopulationLimit`.
For example: let's add these two parameters to the input file above::

  windowType = "phononTransport"
  windowPopulationLimit = 1.0e-6
  temperatures = [3.]
  qMesh = [40,40,40]

Here, we are discarding all phonon states whose equilibrium occupation number is smaller than 1.0e-6.
This will therefore discard, for example, phonon modes away from the Gamma point or optical modes that are too high in energy to be thermally excited at low temperatures.
The resulting calculation will be much faster.
As a result, we can increase the values of `qMesh`, so that we can accurately sample the points close to the Gamma point.




Parallelization
---------------

For this calculation, the bottleneck is typically the construction of the scattering matrix (or the evaluation of a scattering matrix-vector product).
We have three different parallelization schemes.

If you are not familiar with parallelization techniques, check (at least) this https://en.wikipedia.org/wiki/OpenMP and that https://en.wikipedia.org/wiki/Message_Passing_Interface .

* MPI parallelization. We distinguish two cases.
  If we want to compute the action of matrix :math:`\sum_{k'b'} A_{k,k',b,b'} f_{k'b'}`, we MPI-distribute over rows of wavevectors to achieve the best performance. If we want to store the matrix in memory, we parallelize over pairs of wavevectors, using the SCALAPACK layout.

* The calculation of the phonon-phonon coupling is accelerated with Kokkos. Depending on your architecture and installation parameters, Kokkos' code will either run on GPU or on the CPU with OpenMP acceleration. In the former case, remember to set the environmental variable `"export MAXMEM=4"` in the job submission script, or in the command line, to set the available GPU on-board memory (4Gb in this example).

* The summations over band indices when computing the scattering rates is accelerated using OpenMP.

A short guideline to optimize parameters:

* Set the number of MPI processes equal to the number of computing nodes you are requesting.

* Set the number of OpenMP threads equal to the number of physical cores available on each computing node.

* Compile phoebe with Kokkos if you have a GPU. If you do so, make sure that the number of GPU you are using matches the number of MPI processes.


