Running Phoebe
==============

Command line options
--------------------

The basic syntax for launching Phoebe is::

  ./path_to/phoebe -in inputFile -out outputFile

where inputFile.in is the name of the input file with the user-defined variables, and outputFile.out is the name of the output file on which the standard output will be written.


Parallel runs
-------------

Phoebe uses MPI, OpenMP and GPU parallelizations schemes.
It is recommended to compile and run the code with MPI and OpenMP, additionally with CUDA if GPUs are available.
The code can be launched as::

  export OMP_NUM_THREADS=NUM_OMP
  mpirun -np NUM_MPI --bind-to none ./path_to/phoebe -in inputFile -out outputFile

where you should place an appropriate number for NUM_OMP and NUM_MPI
As a first guess (before you optimize it), we recommend to set the number of OpenMP threads equal to the number of cores on a single node, and we recommend to set the number of MPI processes equal to the number of nodes you want to use.

Example:
Suppose we want to run Phoebe on 2 nodes, where each of these nodes comes with 8 cores.
Then, a temptative guess for the parallelization setup is::

  export OMP_NUM_THREADS=8
  mpirun -np 2 --bind-to none ./path_to/phoebe -in inputFile -out outputFile

In this example, Phoebe will be launched with 2 MPI processes, and 8 OpenMP threads for each MPI process.

.. note::
  the 'best' parallelization setup may differ from this first guess we are suggesting, and it is highly dependent on the details of your machine.


GPU acceleration
----------------

Besides making sure that the code is compiled with CUDA, the command line instructions are very similar to the previous case.
An important additional variable to use is `MAXMEM`, which describes the amount of memory available on the GPU in Gigabytes.
The default is set to 16Gb.

*Example*: 
Suppose we want to run Phoebe on 2 nodes, where each of these nodes comes with 8 cores and 1 GPU with 4Gb of GPU on-board memory.
A temptative guess for the parallelization setup is::

  export OMP_NUM_THREADS=8
  export MAXMEM=4
  mpirun -np 2 --bind-to none ./path_to/phoebe -in inputFile.in -out outputFile.out

The variable `MAXMEM` is used to store as many results as possible on the GPU; when the whole GPU memory is filled we return partial results to the CPUs for further processing.
Small test examples should not be impacted by this variable.
See Interaction3Ph for more details.
