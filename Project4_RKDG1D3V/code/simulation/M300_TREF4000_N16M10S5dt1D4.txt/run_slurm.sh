#!/bin/bash
#SBATCH -J dgv1d3vmpi0
#SBATCH -o dgv1d3vmpi0.o%j
#SBATCH -e dgv1d3vmpi0.e%j
#SBATCH -N 1
#SBATCH --ntasks-per-node 28
#SBATCH -p RM
#SBATCH -t 3:00:00
#SBATCH -A ms560hp 
#SBATCH --mail-user=alexander.alekseenko@csun.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

export OMP_NUM_THREADS=1

export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

cd /pylon5/ms560hp/alekseen/DGV0/

mpirun -np $SLURM_NTASKS ./rkdg1d3vMPI.a  >M3ESN16k4M10s5dt1D4.txt 

