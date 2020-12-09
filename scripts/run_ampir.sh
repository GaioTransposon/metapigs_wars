#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=24:00:00
#PBS -l mem=40g
#PBS -N run_ampir
#PBS -M daniela.gaio@student.uts.edu.au

source activate ampir_env

PATH="/shared/homes/12705859/miniconda3/envs/ampir_env/bin:$PATH"

# important to direct to the right R 
R < /shared/homes/12705859/ampir/run_ampir.R --no-save 








