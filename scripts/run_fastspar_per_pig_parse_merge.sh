#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=24:00:00
#PBS -l mem=40g
#PBS -N run_fastspar_per_pig_parse_merge
#PBS -M daniela.gaio@student.uts.edu.au

source activate R36

PATH="/shared/homes/12705859/miniconda3/envs/ampir_env/bin:$PATH"

# important to direct to the right R 
R < /shared/homes/s1/pig_microbiome/fastspar/fastspar_in_per_pig/fastspar_per_pig_parse_merge.R --no-save 
