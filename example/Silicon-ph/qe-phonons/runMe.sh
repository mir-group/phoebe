#!/usr/bin/env bash

export QE_PATH="/your/path/to/QE/bin"
export NMPI=4
export NPOOL=4

mpirun -np $NMPI $QE_PATH/pw.x -npool $NPOOL -in scf.in > scf.out

mpirun -np $NMPI $QE_PATH/ph.x -npool $NPOOL -in ph.in > ph.out

$QE_PATH/q2r.x -in q2r.in > q2r.out
