#!/bin/bash
#PBS -l ncpus=20
#PBS -l walltime=120:00:00
#PBS -l mem=40g
#PBS -N run_ampir
#PBS -M daniela.gaio@student.uts.edu.au

source activate ampir_env

export PATH=/shared/homes/12705859/miniconda3/condabin:/shared/homes/s1/pig_microbiome/mmseqs_nextflow/mmseqs2/bin:/shared/homes/s1/pig_microbiome/mmseqs_nextflow:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/pbs/bin:/shared/homes/12705859/git/bwa:/shared/homes/12705859/drep/bin:/shared/homes/12705859/edirect:/shared/homes/12705859/bin:/shared/homes/12705859/.aspera/connect/bin

# important to direct to the right R 
/shared/homes/12705859/miniconda3/envs/ampir_env/bin/R < /shared/homes/12705859/ampir/run_ampir.R --no-save 








