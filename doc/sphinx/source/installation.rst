Installing Phoebe
=================

Download
--------

The code is available (free of charge with an open-source MIT license) at its [github page](https://github.com/mir-group/phoebe).
If checking out from the GitHub repository, make sure to use the master branch for production. The developer's branch, while we will attempt to keep it in working condition, is not recommended for public usage.


Prerequisites
-------------

The installation requires the following packages pre-installed:

* CMake (a recent version);
  
* a C++ compiler with C++14 support. We regularly test with GCC, but Intel and Clang should work too;
  
* MPI (although the code can compile without);
  
* Optional: OpenMP;
  
* Optional: CUDA (for GPU acceleration);
  
* Internet connection, to download external libraries.



Cmake build
-----------

Basic build
^^^^^^^^^^^

To install Phoebe, type::

  git submodule update --init
  mkdir build
  cd build
  cmake ..
  make -j$(nproc)

where you should substitute `nproc` with a number of parallel compilation jobs.
This will create the executable `phoebe` in the `build` directory.

CMake will inspect the paths found in the environmental variable `LD_LIBRARY_PATH` to verify the existence of an installed copy of the SCALAPACK library and link it. If not found, the installation will compile a copy of the SCALAPACK.

HDF5 build
^^^^^^^^^^

Phoebe can make use of HDF5 through the HighFive library to write the electron-phonon matrix elements in the elPhQeToPhoebe app, 
as well as any app which reads in and uses these matrix elements. 
This is highly recommended, as it speeds up what can be time consuming IO, and also significantly reduces file sizes. 
When built using cmake with the flag `-DHDF5_AVAIL=ON`, Phoebe will be built with HDF5. If MPI is also present, 
Phoebe will be built to perform HDF5 operations in parallel. 

If, for some reason, a user has MPI present, but has built a copy of HDF5 which does not link to MPI and therefore cannot 
perform parallel read/write operations, they must build Phoebe using the -DHDF5_SERIAL=ON cmake option to force serial HDF5 operations.

Note: When building on Ubuntu, one may need to specify the location of hdf5 to cmake. This can be done, for example, when using 
libhdf5-openmpi-dev::

  cmake .. -DCMAKE_CXX_STANDARD_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/" -DCMAKE_CXX_FLAGS="-I/usr/include/hdf5/openmpi/"


OpenMP build
^^^^^^^^^^^^

::

  git submodule update --init
  mkdir build
  cd build
  cmake .. -DKokkos_ENABLE_OPENMP=ON -DOMP_AVAIL=ON
  make -j$(nproc)


OpenMP + CUDA build
^^^^^^^^^^^^^^^^^^^

::
   
  git submodule update --init
  mkdir build
  cd build
  cmake .. -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_VOLTA70=ON -DOMP_AVAIL=ON
  make -j$(nproc)

Replace VOLTA70 (e.g. for V100s GPUs) with the arch of your GPU.

Note on Kokkos
^^^^^^^^^^^^^^
Phoebe utilizes Kokkos to generate either OpenMP or CUDA code to accelerate parts of the code.
Thus Phoebe accepts all the CMake arguments of Kokkos, which can improve performance.
For example, you should specify -DKokkos_ARCH_KNL=ON when building for Knight's Landing nodes.





Compiling the documentation
---------------------------

In order to compile the documentation, you need to have installed on your machine:

* doxygen
  
* graphviz
  
* pdflatex (to render equations)
  
Typically for Unix machines, these packages are commonly found on package managers (apt, pacman, brew, ...).

Then type::

  cd build
  make doc

Note that compiling the documentation doesn't require compiling the code.



Installation instructions for common workstations and supercomputers
--------------------------------------------------------------------

Ubuntu
^^^^^^

To install (without GPU support)::

  sudo apt install cmake gcc doxygen graphviz libomp-dev libopenmpi3 libhdf5-openmpi-dev
  git submodule update --init
  mkdir build
  cd build
  cmake .. -DKokkos_ENABLE_OPENMP=ON -DOMP_AVAIL=ON -DCMAKE_CXX_STANDARD_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/" -DCMAKE_CXX_FLAGS="-I/usr/include/hdf5/openmpi/"`
  make -j$(nproc)
  make doc

Note that paths to the hdf5 library may need to be updated
Tested on Ubuntu 20.04.

MacOs
^^^^^

We have experienced troubles linking the SCALAPACK library, especially when linking it together with the libgfortran library.
If libgfortran is not found, try adding it specifically to LD_LIBRARY_PATH::

  export LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/libgfortran/

In particular, if you are using a version of gcc installed using homebrew, you might need to link the `Cellar` copy of libgfortran. As an example working for gcc v9.3.0_1 is::

  export LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/Cellar/gcc/9.3.0_1/lib/gcc/9/) 

