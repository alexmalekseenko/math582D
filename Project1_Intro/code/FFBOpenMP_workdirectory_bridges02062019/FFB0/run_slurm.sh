#!/bin/bash
#SBATCH -J ffb0
#SBATCH -o "ffb0.%j.%N.out"
#SBATCH -p RM
#SBATCH -N 1
#SBATCH --ntasks-per-node 28
#SBATCH -t 7:00:00
#SBATCH -A mp5fp1p 
#SBATCH --mail-user=alexander.alekseenko@csun.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

export OMP_NUM_THREADS=28

cd /pylon5/mp5fp1p/alekseen/FFB0/

./ffbM650.a >outM650ffbM41tr3_TREF10000_Decomp_dt2D6.txt
