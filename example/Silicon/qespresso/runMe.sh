#!/usr/bin/env bash

pw.x -in scf.in > scf.out

ph.x -in ph.in > ph.out

q2r.x -in q2r.in > q2r.out

pw.x -in scf.in > scf.out

pw.x -in nscf.in > nscf.out

wannier90.x -pp si

pw2wannier90.x -in pw2wan.in > pw2wan.out

wannier90.x si
