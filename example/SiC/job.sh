#!/bin/bash
#   *** Single Serial Job in v100 Queue ***
#
# Notes:
#
#   -- Copy/edit this script as desired.  Launch by executing
#      "sbatch sample.slurm" on a Longhorn login node.
#
#   -- Serial codes run on a single node (upper case N = 1).
#        A serial code ignores the value of lower case n,
#        but slurm needs a plausible value to schedule the job.
#----------------------------------------------------
#SBATCH -J myjob           # Job name
#SBATCH -p v100            # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH -A DMR20009       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
pwd
date
module load gcc/7.3.0 cuda/10.2
module list

# Launch serial code...
CUDA_LAUNCH_BLOCKING=1 ../../cudabuild/phoebe -in input.in
