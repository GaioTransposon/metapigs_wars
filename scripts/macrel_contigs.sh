
## usage example: 
## export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_141.fa
## export out_dir=/shared/homes/12705859/macrel_work/out_contigs141
## qsub -V macrel_contigs.sh


#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=72:00:00
#PBS -l mem=40g
#PBS -N macrel_contigs
#PBS -M daniela.gaio@student.uts.edu.au

source activate macrel_env

cd /shared/homes/12705859/macrel_work

macrel contigs \
    --fasta $input_file \
    --output $out_dir \
    -t 8