#!/bin/bash

#PBS -P vf71
#PBS -q normalbw
#PBS -l ncpus=84
#PBS -l mem=756GB
#PBS -l jobfs=200GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/vf71+gdata/hh5
#PBS -M lilian.fierroarcos@utas.edu.au
#PBS -m abe
#PBS -l wd

module use /g/data/hh5/public/modules
module load conda/analysis3
python3 run_dbpm_gridded.py
