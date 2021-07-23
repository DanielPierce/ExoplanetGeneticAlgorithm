#!/bin/bash
#SBATCH --job-name="Exoplanet GA Testing"
#SBATCH --output="testdata.%j.%N.out"
#SBATCH --error="testdata.%j.%N.err"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --export=ALL
#SBATCH --account=fsc102

## Setup run environment
module reset
module load gcc/10.2.0
module load python

## Run script (replace with real)
python test_exp.py