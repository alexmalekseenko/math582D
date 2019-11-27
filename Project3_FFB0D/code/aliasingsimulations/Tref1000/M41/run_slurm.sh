#!/bin/bash
#SBATCH -J ffb1
#SBATCH -o "ffb1.%j.%N.out"
#SBATCH -p RM
#SBATCH -N 1
#SBATCH --ntasks-per-node 28
#SBATCH -t 00:20:00
#SBATCH -A ms560hp 
#SBATCH --mail-user=alexander.alekseenko@csun.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

export OMP_NUM_THREADS=28

cd /pylon5/ms560hp/alekseen/FFB1/

./M3alisng.a >outtestM3aliasingTref1000collisionop_M41.txt
