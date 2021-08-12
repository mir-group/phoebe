Running Phoebe
==============

Command line options
--------------------

The basic syntax for launching Phoebe is::

  ./path_to/phoebe -in inputFile -out outputFile

where ``inputFile.in`` is the name of the input file with the user-defined variables, and ``outputFile.out`` is the name of the output file to which output will be written.


Parallel runs
-------------

Phoebe uses MPI, OpenMP and GPU parallelizations schemes.
It is recommended to compile and run the code with MPI and OpenMP, additionally with CUDA (through Kokkos) if GPUs are available.
The code can be launched as::

  export OMP_NUM_THREADS=NUM_OMP
  mpirun -np NUM_MPI --bind-to none ./path_to/phoebe -in inputFile -out outputFile

where you should place an appropriate number for NUM_OMP and NUM_MPI.
As a first attempt (before you optimize it), we recommend to set the number of OpenMP threads equal to the number of cores on a single node, and we recommend to set the number of MPI processes equal to the number of nodes you want to use. This will maximize the memory available to each MPI process, as the OMP threads on each node will share the entirety of each node's memory.

**Example:**
Suppose we want to run Phoebe on 2 nodes, where each of these nodes comes with 8 cores.
Then a first guess for the parallelization setup is::

  export OMP_NUM_THREADS=8
  mpirun -np 2 --bind-to none ./path_to/phoebe -in inputFile -out outputFile

In this example, Phoebe will be launched with 2 MPI processes, and 8 OpenMP threads for each MPI process. This information is printed at the start of each Phoebe run. Often, other good configurations can be found with more MPI processes per node -- just be sure that the memory shared among the remaining OMP threads/per core is enough to run the calculation.

.. note::
  The best parallelization setup will likely differ from this suggestion, as it is highly dependent on the details of your machine.
  
GPU acceleration
----------------

To run the parts of Phoebe which are optimized for use on GPUs (the calculations of phonon-phonon and electron-phonon scattering rates), first ensure that the code is compiled with CUDA through Kokkos, as described on the :ref:`installation` page. The command line instructions are very similar to the previous case, barring the inclusion of one important additional variable, ``MAXMEM``, which describes the amount of memory available on the GPU in gigabytes.
The default is set to 16GB.

**Example**:
Suppose we want to run Phoebe on 2 nodes, where each of these nodes comes with 8 cores and 1 GPU with 4GB on-board memory.
A tentative guess for the parallelization setup is::

  export OMP_NUM_THREADS=8
  export MAXMEM=4
  mpirun -np 2 --bind-to none ./path_to/phoebe -in inputFile.in -out outputFile.out

The variable ``MAXMEM`` is used to store as many results as possible on the GPU; when the whole GPU memory is filled, Phoebe returns partial results to the CPUs for further processing. Small test examples should not be impacted by this variable.
See the Phoebe run in the :ref:`phononTransport` for more details.

.. note::
   The number of MPI processes must be equal to the number of GPU being used.

  
Electron-phonon parallelization
-------------------------------

In order to perform the Wannier interpolation of the electron-phonon coupling, the coupling tensor describing this interaction muts fit in the memory available to each single MPI process and, if applicable, in the GPU memory.
However, the electron-phonon tensor may be extremely large (see the size of the *.phoebe.elph.hdf5 or *.phoebe.elph.dat files for a guess on the size), causing out-of-memory errors.

To overcome this, the tensor can be distributed over different MPI processes, reducing the memory requirements of the calculation in exchange for a slowdown due to increased communications between MPI nodes.
The parallelization of this tensor is performed over the :math:`\boldsymbol{R}_e` index (the conjugate variable of the k-point index, see the theory section).
This parallelization is controlled by the `poolSize` or `ps` command line parameter, where poolSize is the number of MPI processes over which the coupling tensor is distributed.

For example, suppose we want to run Phoebe on a computer with four nodes, each with eight-core CPUs and one GPU.
As mentioned above, the number of MPI processes must be equal to the number of GPU available, i.e. 4.
Next, we want to use all CPUs present on a node, therefore we set OMP_NUM_THREADS to 8.
Let's suppose that the electron-phonon coupling tensor fits in the combined memory of 2 GPUs.
Therefore, we distribute the coupling tensor over a pool of 2 MPI processes, as::
  
  export OMP_NUM_THREADS=8
  mpirun -np 2 --bind-to none ./path_to/phoebe --poolSize 2 -in inputFile.in -out outputFile.out

Or equivalently::
  
  export OMP_NUM_THREADS=8
  mpirun -np 2 --bind-to none ./path_to/phoebe -ps 2 -in inputFile.in -out outputFile.out

.. note:: the poolSize parameter must be an integer divisor of the number of MPI processes.
For example, if using 10 MPI processes, the poolSize can only be set to 1, 2, 5, or 10, with the highest number having the lowest memory footprint but slowest performance.
