#!/usr/bin/env bash

# Note: QE executables come from the patched version of QE
# Wannier90 is unmodified

mpirun -np 4 pw.x -npool 4 -in scf.in > scf.out
mpirun -np 4 ph.x -npool 4 -in ph.in > ph.out
q2r.x -in q2r.in > q2r.out

#mpirun -np 1 pw.x -npool 1 -in scf.in > scf.out
#mpirun -np 1 pw.x -npool 1 -in nscf.in > nscf.out
#~/software/wannier90-3.0.0/wannier90.x -pp si
#mpirun -np 1 pw2wannier90.x -in pw2wan.in > pw2wan.out
#~/software/wannier90-3.0.0/wannier90.x si

#mpirun -np 4 pw.x -npool 4 -in scf.in > scf.out
#mpirun -np 4 pw.x -npool 4 -in nscf.in > nscf.out
#~/software/wannier90-3.0.0/wannier90.x -pp si
#mpirun -np 4 pw2wannier90.x -in pw2wan.in > pw2wan.out
#~/software/wannier90-3.0.0/wannier90.x si

#mpirun -np 4 ~/git/q-e/bin/pw.x -npool 4 -in scf.in > scf.out
#mpirun -np 4 ~/git/q-e/bin/pw.x -npool 4 -in nscf.in > nscf.out
#~/software/wannier90-3.0.0/wannier90.x -pp si
#mpirun -np 4 ~/git/q-e/bin/pw2wannier90.x -in pw2wan.in > pw2wan.out
#~/software/wannier90-3.0.0/wannier90.x si
