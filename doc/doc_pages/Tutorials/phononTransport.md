@page phTrTut Phonon Transport Tutorial

@section Synopsis Synopsis
<ol>
<li> Run QE `pw.x` as an `scf` calculation </li>
<li> Run QE `ph.x` phonon calculation </li>
<li> Run QE `q2r.x` to process the dynamical matrix </li>
<li> Compute 3rd order energy derivatives </li>
<li> Run Phoebe's app @ref phtr </li>
</ol>

In this tutorial, we want to compute the lattice thermal conductivity of Silicon.
We will use Quantum ESPRESSO as the code that provides DFT parameters.
As described in the @ref Theory section of this manual, we need Quantum ESPRESSO to compute the phonon dynamical matrix and the third derivative of the total energy with respect to ionic displacements.

Note that we assume the reader to be familiar with Quantum ESPRESSO.
Several tutorials can be found on <a href="https://www.quantum-espresso.org/resources/tutorials">Quantum ESPRESSO's website</a>, which cover more complicated DFT calculations than that described here.


@section step1 Step 1: Pw
First, we need to compute the total energy of the silicon crystal unit cell.
This calculation will create the ground state charge density and wavefunctions that are needed for later.

To run, go to the folder `./example/Silicon/qespresso` in the phoebe repository.
The file `scf.in` is the input file for the `pw.x` executable.
The input for a total energy DFT calculation of Quantum ESPRESSO for a silicon crystal, is
~~~~~~~~~~{.c}
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
~~~~~~~~~~
A detailed description of all this parameters can be found on <a href="https://www.quantum-espresso.org/Doc/INPUT_PW.html">Quantum ESPRESSO's website</a>.
The most important parameters to be tweaked and modified in a research project are
<ul>
<li> `K_POINTS`: parameter controlling the integration mesh of wavevectors in the Brillouin zone. Phonon properties should be converged against this mesh (more wavevectors is better). Tip: `ph.x` calculations are faster when the k-mesh is gamma-centered. </li>
<li> `ecutwfc`: parameter controlling the number of G-vectors used in the plane-wave expansion of the wavefunction. Phonon frequencies should be checked against this value. </li>
<li> `conv_thr` </li> this parameter controls total energy convergence. Note that a poorly converged conv_thr may result in poorly converged phonon properties.
<li> `prefix`: prefix of some output files. Make sure to use a consistent value of prefix throughout your calculations.
<li> `outdir`: name of the scratch folder. Must be used consistently throughout thre run, so to point to the correct files. </li>
</ul>
This list is obviously not complete, and for your research project you may need to use more functionalities of QE's `pw.x`.

Simply run it as

    /path/to/qe/pw.x -in scf.in > scf.out

after substituting the suitable path to the `pw.x` executable.


- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


@section step2 Step 2: Ph
The input file `ph.in` is as follows:
~~~~~~{.c}
phonons of Si
 &inputph
  tr2_ph = 1.0d-14,
  prefix = 'silicon',
  ldisp = .true.,
  nq1=4, nq2=4, nq3=4
  outdir = "./out",
  fildyn = 'silicon.dyn',
 /
~~~~~~
The values of `nqX` select the Monkhorst-Pack grid of q-points centered at Gamma, for which we will compute the phonon properties.
Here it's important that `prefix` and `outdir` are the same as those used in the `pw.x` calculation of before.
Use a good value of `tr2_ph` (smaller is better, but harder to converge), which (indirectly) checks the convergence of phonon frequencies.

Run the code as

    /path/to/qe/bin/ph.x -in ph.in > ph.out

If the code executes correctly and completely, you should see a number of files called `{fildyn}*`, as many files as the number of irreducible q-points.

<blockquote>
We strongly recommend to parallelize this calculation, and to use `npool`. This parameter controls the parallelization over k-vector. It gives an almost linear speedup, but is bound by the number of k-points. For example:
~~~~~~~{.c}
      mpirun -np 4 /path/to/qe/bin/ph.x -npool 4 -in ph.in > ph.out
~~~~~~~
</blockquote>


<blockquote>
This is a simple example built such that the parameters used here yield reasonable results.
In general, we strongly recommend to test convergence of the phonon frequencies with respect to the k-point mesh, the q-point mesh, the wavefunction cutoff and the `tr2_ph` parameter.
We also recommend to use a small `conv_thr`.
</blockquote>


- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


@section step3 Step 3: Q2r

The code ph.x has created the `silicon.dyn*` files, which contain the dynamical matrix at every irreducible q-point.
Now, we run `q2r.x` to Fourier transform the dynamical matrices in the reciprocal space representation to the real space representation, where they represent the interatomic force constants.
The input file `q2r.in` is minimal:
~~~~~~{.c}
 &input
   fildyn='silicon.dyn',
   flfrc='silicon.fc'
 /
~~~~~~
where the first variable must match the path to the dynamical matrices set earlier in `ph.x`, and `flfrc` is the output file with the force constants. 

In the working folder `./example/Silicon/qespresso` run the command

    ./path/to/qe/bin/q2r.x -in q2r.in > q2r.out

If the code run successfully, you should see a new file `silicon.fc`.


- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


@section step4 Step 4: 3rd derivatives
Blank!


- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


@section step5 Step 5: Phoebe

The typical input file looks like this:
~~~~~~~~~~~{.c}
appName = "phononTransport"
phD2FileName = "./qespresso/silicon.fc",
phD3FileName = "./qespresso/silicon.fc3"
sumRuleD2 = "crystal"
qMesh = [10,10,10]
temperatures = [300.]
smearingMethod = "adaptiveGaussian"
solverBTE = ["variational"]
~~~~~~~~~~~

Let's go through this parameters one by one:

<ol>
<li> @ref appName = `"phononTransport"` triggers the calculation of phonon transport properties </li> 
<li> @ref phD2FileName must point to the `flfrc` file produced by `q2r.x` </li> 
<li> @ref phD2FileName must point to the file of third derivatives  </li>
<li> @ref sumRuleD2 allows us to re-enforce the translational-invariance of force constants, that is broken by numerical errors. After imposing this sum rule, acoustic phonon frequencies to go to zero at the gamma point. </li>
<li> @ref qMesh is the size of the grid of wavevectors used to integrate the Brillouin zone. Note that the value used here is very unconverged, so that the example can finish in a short amount of time.
<blockquote>
Results must be converged against values of `qMesh`!
</blockquote>
</li>
<li> @ref temperatures sets the temperature in Kelvin  </li>
<li> @ref smearingMethod sets the algorithm to approximate the Dirac-delta conserving energy. Using the "adaptiveGaussian" scheme is particular convenient as the width of the gaussian is automatically adjusted. With "gaussian" scheme instead, you should converge the @ref smearingWidth parameter together with the @ref qMesh.   </li> 
<li> @ref solverBTE Selects the algorithm to solve the linearized Boltzmann equation. If not specified, we only compute results within the relaxation time approximation. Here, we are using the variational solver to find the solution to the BTE. </li> 
</ol>

With this input, we can compute the phonon contribution to thermal conductivity of silicon.

<blockquote>
By default, isotopic scattering at natural abundances is included in the scattering matrix. To disable or modify it, check the parameters @withIsotopeScattering and @massVariance.
</blockquote>

<blockquote>
In several studies you may want to include boundary scattering. To include it, use the parameter @ref boundaryLength.
</blockquote>




@section highmem Tradeoffs between speed and memory
There's a parameter @ref scatteringMatrixInMemory that you need to consider.
If we set this parameter to true, we store the scattering matrix in memory.
If false, we only compute the action of the scattering matrix, without ever storing all of it in memory.

There is no `best` choice here, rather, you should decide what's best for your case and decide which tradeoff works best for you.

<ul>
<li>
Option 1: @ref scatteringMatrixInMemory = true. The scattering matrix occupies \f$ 16 (3 N_{atoms} N_{q-points})^2 / 1024^3 \f$ Gygabytes, if no window is used. This number can be pretty large (even Terabytes), and you should make sure that your HPC allocation has enough memory for storing this large matrix. Given the size, we only allow you to run the code with a single temperature.

In exchange, iterative or variational solvers of the BTE are extremely cheap, and the cost of your simulation is basically the cost of constructing the scattering matrix. Moreover, you get access to the @ref solverBTE="relaxons" type of BTE solver.
</li>

<li>
Option 2: @ref scatteringMatrixInMemory = false. The memory footprint is much lighter (the square root of before), so that the same calculation can be run on fewer CPUs. You can compute the thermal conductivity for multiple temperatures in the same run. The calculation of properties within the relaxation time approximation is as expensive as above (if this is what you care about, definitely use less memory).

In exchange, iterative or variational BTE solvers are much slower.
In fact, at each iteration you need to recompute the scattering matrix.
The cost of the calculation therefore grows linearly with the number of iterations of the iterative solver (which may be significant).
You also cannot diagonalize the scattering matrix with @ref solverBTE="relaxons".
</li>
</ul>





@section lowT Low temperature thermal conductivity
At low temperatures, only phonons with small energies are thermally excited and most states are empty.
However, if we use the same input file as the one above, we are sampling all phonon states.
As a result, we end up spending a lot of time computing phonon states that don't contribute to transport.

To address this, we have the parameters @ref windowType, @ref windowEnergyLimit, and @windowPopulationLimit.
For example: let's add these two parameters to the input file above:
~~~~~~~~~~~{.c}
windowType = "phononTransport"
windowPopulationLimit = 1.0e-6
temperatures = [3.]
qMesh = [40,40,40]
~~~~~~~~~~~
Here, we are discarding all phonon states whose equilibrium occupation number is smaller than 1.0e-6.
This will therefore discard, for example, phonon modes away from the Gamma point or optical modes that are too high in energy to be thermally excited at low temperatures.
The resulting calculation will be much faster.
As a result, we can increase the values of `qMesh`, so that we can accurately sample the points close to the Gamma point.




@section para Parallelization
For this calculation, the bottleneck is typically the construction of the scattering matrix (or the evaluation of a scattering matrix-vector product).
We have three different parallelization schemes.

If you are not familiar with parallelization techniques, check (at least) <a href="https://en.wikipedia.org/wiki/OpenMP"> this</a> and <a href="https://en.wikipedia.org/wiki/Message_Passing_Interface"> this</a> to begin with.

<ol>
<li>
MPI parallelization. We distinguish two cases.
If we want to compute the action of matrix \f$ \sum_{k'b'} A_{k,k',b,b'} f_{k'b'} \f$, we MPI-distribute over rows of wavevectors to achieve the best performance. If we want to store the matrix in memory, we parallelize over pairs of wavevectors, using the SCALAPACK layout.
</li>
<li>
The calculation of the phonon-phonon coupling is accelerated with Kokkos. Depending on your architecture and installation parameters, Kokkos' code will either run on GPU or on the CPU with OpenMP acceleration. In the former case, remember to set the environmental variable `"export MAXMEM=4"` in the job submission script, or in the command line, to set the available GPU on-board memory (4Gb in this example).
</li>
<li>
The summations over band indices when computing the scattering rates is accelerated using OpenMP.
</li>
</ol>

A short guideline to optimize parameters:
<ul>
<li>
Set the number of MPI processes equal to the number of computing nodes you are requesting.
</li>
<li>
Set the number of OpenMP threads equal to the number of physical cores available on each computing node.
</li>
<li>
Compile phoebe with Kokkos if you have a GPU. If you do so, make sure that the number of GPU you are using matches the number of MPI processes.
</li>
</ul>