export QE_PATH=/path/to/phoebe-patched-qe/bin/
export PHOEBE_PATH=/path/to/phoebe/build/
# increase these if possible to speed up the calculation
export NMPI=4
export NPOOL=4

# run the phonon calculation
mpirun -np $NMPI $QE_PATH/pw.x -npool $NPOOL -in scf.in > scf.out
mpirun -np $NMPI $QE_PATH/ph.x -npool $NPOOL -in ph.in > ph.out
$QE_PATH/q2r.x -in q2r.in > q2r.out

# run the wannier calculation
mpirun -np $NMPI $QE_PATH/pw.x -npool $NPOOL -in nscf.in > nscf.out
$QE_PATH/wannier90.x -pp cu
$QE_PATH/pw2wannier90.x -in pw2wan.in > pw2wan.out
mpirun -np $NMPI $QE_PATH/wannier90.x cu

# run phoebe to convert qe -> phoebe.hdf5 file
mpirun -np $NMPI $PHOEBE_PATH/phoebe -in qeToPhoebeWannier.in > 2phoebe.out
