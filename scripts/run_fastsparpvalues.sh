
## usage: 
## export your_dir=/shared/homes/s1/pig_microbiome/fastspar
## export input_file=fastsparGTDB_in.txt
## qsub -V run_fastsparpvalues.sh


#!/bin/bash
#PBS -l ncpus=12
#PBS -l walltime=72:00:00
#PBS -l mem=10g
#PBS -N run_fastsparpvalues
#PBS -M daniela.gaio@student.uts.edu.au

source activate fastspar_env

cd $your_dir

# STEP 4: fastspar pvalues (pvalues are calculated from the correlations)
fastspar_pvalues --threads 10 --otu_table $input_file --correlation median_correlation.tsv --prefix bootstrap_correlation/cor_ --permutations 1000 --outfile pvalues.tsv
