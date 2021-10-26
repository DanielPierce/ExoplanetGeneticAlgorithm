#!/bin/bash
#SBATCH --job-name="quick test"
#SBATCH --output="testdata.%j.%N.out"
#SBATCH --error="testdata.%j.%N.err"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --export=ALL
#SBATCH --account=TG-CCR180031
#SBATCH -t 00:10:00
#SBATCH --mail-type=BEGIN,END,FAIL

## Setup run environment
module reset
module load gcc/10.2.0
module load python
module load py-pip
module load py-matplotlib
module load py-numpy

## Run script
pip install --user python-dateutil
pip install --user lightkurve
pip install --user astropy
pip install --user deap
pip install --user jsonpickle
pip install --user mpmath
pip install --user orbitalpy
pip install --user openpyxl

echo --+--PIP LIST--+--
pip list

python -m PSEUGA
