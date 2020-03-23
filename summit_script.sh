#!/bin/bash
#BSUB -P MY_ACCT_NUM
#BSUB -W 0:20
#BSUB -nnodes 1
#BSUB -J INCFLO
#BSUB -o INCFLOo.%J
#BSUB -e INCFLOe.%J
module load gcc
module load cuda/9.1.85
module list
set -x
omp=1
export OMP_NUM_THREADS=${omp}
EXE="./incflo3d.gnu.MPI.CUDA.ex"
jsrun -n 1 -a 1 -g 1 -c 1 --bind=packed:${omp} ${EXE} inputs > output.txt
