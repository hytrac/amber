#!/bin/bash
#SBATCH --job-name=amber
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=240GB
#SBATCH -t 24:00:00
#SBATCH -o std.log

# Set environment variables
export KMP_STACKSIZE=128m

# Go to the job scratch directory
cd $SLURM_SUBMIT_DIR

# Compile and run job
make clean; make
./amber.x < input/input.txt > output/log
