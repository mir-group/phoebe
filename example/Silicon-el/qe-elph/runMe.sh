#!/usr/bin/env bash

export QE_PATH=/your/path/to/patched/phoebe-quantum-espresso/bin
export NMPI=4
export NPOOL=4

# Note: must point to the patched version of QE
# (you can use an unmodified version of Wannier90)

mpirun -np $NMPI $QE_PATH/pw.x -npool $NPOOL -in scf.in > scf.out
mpirun -np $NMPI $QE_PATH/ph.x -npool $NPOOL -in ph.in > ph.out
$QE_PATH/q2r.x -in q2r.in > q2r.out

mpirun -np $NMPI $QE_PATH/pw.x -npool $NPOOL -in nscf.in > nscf.out
$QE_PATH/wannier90.x -pp si
mpirun -np $NMPI $QE_PATH/pw2wannier90.x -in pw2wan.in > pw2wan.out
$QE_PATH/wannier90.x si
