Installing Phoebe
=================

Download
--------

The code is available (free of charge with an open-source MIT license) at its `github page <https://github.com/mir-group/phoebe>`__.
When checking out from the GitHub repository, make sure to use the master branch. The other branches may be functional, but are not guaranteed to work and should be used with caution.

Prerequisites
-------------

The installation requires the following packages are available:

* CMake (a recent version)

* A C++ compiler with C++14 support, whether GCC, Intel, or Clang

* MPI (although the code can compile without)

* Optional: OpenMP

* Optional: CUDA (for GPU acceleration)

* An internet connection, to download external libraries


Cmake build
-----------

Basic build
^^^^^^^^^^^

To install Phoebe, type::

  git submodule update --init
  mkdir build
  cd build
  cmake ..
  make -j $(nproc)

where you should substitute ``nproc`` with the number of cores available for parallel compilation. This will create the executable ``phoebe`` in the ``build`` directory.

CMake will inspect the paths found in the environmental variable ``LD_LIBRARY_PATH`` to verify the existence of an installed copy of the ScaLAPACK library and link it. If not found, the installation will download and compile a copy of ScaLAPACK.

HDF5 build
^^^^^^^^^^

Phoebe can make use of HDF5 through the `HighFive <https://github.com/BlueBrain/HighFive>`__ library to read and write the electron-phonon matrix elements. This is highly recommended, as it speeds up what can be time consuming I/O, and also significantly reduces file sizes.
When built using CMake with the flag ``-DHDF5_AVAIL=ON``, Phoebe will be built with HDF5. If MPI is also present,
Phoebe will be built to perform HDF5 operations in parallel.

If, for some reason, a user has MPI present, but has built a copy of serial HDF5, which does not link to MPI and therefore cannot
perform parallel read/write operations, they must build Phoebe using the ``-DHDF5_SERIAL=ON`` CMake option to force serial HDF5 operations.

.. note::
  When building on Ubuntu, one may need to specify the location of HDF5 to CMake. This can be done, for example, when using libhdf5-openmpi-dev::

   cmake .. -DCMAKE_CXX_STANDARD_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/" -DCMAKE_CXX_FLAGS="-I/usr/include/hdf5/openmpi/"


OpenMP build
^^^^^^^^^^^^

::

  git submodule update --init
  mkdir build
  cd build
  cmake .. -DKokkos_ENABLE_OPENMP=ON -DOMP_AVAIL=ON
  make -j$(nproc)


Kokkos Build (For GPU use)
^^^^^^^^^^^^^^^^^^^

Replace VOLTA70 (e.g. for V100s GPUs) with the arch of your GPU. Phoebe utilizes Kokkos to generate either OpenMP or CUDA code to accelerate parts of the code. To build for use with Kokkos, specify the flags to enable CUDA and OpenMP functionality.
Phoebe accepts all the CMake arguments of Kokkos, which can improve performance.
For example, you should specify ``-DKokkos_ARCH_KNL=ON`` when building for
Knight's Landing nodes. Additionally, replace VOLTA70 (e.g. for V100s GPUs) as shown below with the arch of your GPU.

::

  git submodule update --init
  mkdir build
  cd build
  cmake .. -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_VOLTA70=ON -DOMP_AVAIL=ON
  make -j$(nproc)

Compiling the documentation
---------------------------

In order to compile the documentation locally (the same documentation as on the Phoebe website), you need to have the following available on your machine:

* doxygen

* graphviz

* pdflatex (to render equations)

Then type::

  cd build
  make doc

Note that compiling the documentation doesn't require compiling the code.


Installation instructions for specific systems
--------------------------------------------------------------------

Ubuntu
^^^^^^

To install (without GPU support)::

  sudo apt install cmake gcc doxygen graphviz libomp-dev libopenmpi3 libhdf5-openmpi-dev
  git submodule update --init
  mkdir build
  cd build
  cmake .. -DKokkos_ENABLE_OPENMP=ON -DOMP_AVAIL=ON -DCMAKE_CXX_STANDARD_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/" -DCMAKE_CXX_FLAGS="-I/usr/include/hdf5/openmpi/"
  make -j$(nproc)
  make doc

Note that paths to the HDF5 library may need to be updated.
Tested on Ubuntu 20.04.

MacOS
^^^^^

* We have encountered difficulty linking the ScaLAPACK library, especially when linking with libgfortran. If libgfortran is not found, try adding it specifically to ``LD_LIBRARY_PATH`` or ``LIBRARY_PATH`` as follows:
  ::

    export LIBRARY_PATH=$LIBRARY_PATH:/path/to/libgfortran/

  In particular, if you are using a version of gcc installed using homebrew, you might need to link the "Cellar" copy of libgfortran. As an example working for gcc v9.3.0_1 is::

    export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/Cellar/gcc/9.3.0_1/lib/gcc/9/)

* Additonally, there exists an issue when building with the Apple Clang compiler
  and the Eigen library, specifically when Eigen is built using OpenMP with a c++ std>11. We recommend either building without OpenMP (``cmake -DOMP_AVAIL=OFF ../``), or using a different compiler.
