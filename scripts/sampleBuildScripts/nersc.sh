
# run this in your build directory
# This worked on nersc for the gcc
# compiler on 10/29/20

module load cmake/3.18.2
module swap PrgEnv-intel PrgEnv-gnu
module load cray-hdf5
export CRAYPE_LINK_TYPE=dynamic
cmake -DOMP_AVAIL=OFF ../  # omp causes problems on nersc currently
make -j 12 

