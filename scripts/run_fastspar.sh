
## usage: 
## export your_dir=/shared/homes/s1/pig_microbiome/fastspar
## export input_file=fastsparGTDB_in.txt
## qsub -V run_fastspar.sh


#!/bin/bash
#PBS -l ncpus=12
#PBS -l walltime=72:00:00
#PBS -l mem=10g
#PBS -N run_fastspar
#PBS -M daniela.gaio@student.uts.edu.au

source activate fastspar_env

cd $your_dir

# STEP 1: setup the bootstrap 
mkdir bootstrap_counts
fastspar_bootstrap --otu_table $input_file --number 1000 --prefix bootstrap_counts/otu_data

# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
mkdir bootstrap_correlation
for item in `ls bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation bootstrap_correlation/cor_$N --covariance bootstrap_correlation/cov_$N
done

# STEP 3: run fastspar
fastspar -y --threads 10 --otu_table $input_file --correlation median_correlation.tsv --covariance median_covariance.tsv --iterations 5 
