#!/bin/bash
#PBS -l ncpus=30
#PBS -l walltime=72:00:00
#PBS -l mem=10g
#PBS -N runFastspEACH
#PBS -M daniela.gaio@student.uts.edu.au

source activate fastspar_env

cd $your_dir

# STEP 1: setup the bootstrap 
for input_file in `ls fastsparGTDB_in_*.txt` 
do
filename=$(basename $input_file)     
N="${filename%.*}"
mkdir $N
mkdir $N/bootstrap_counts
mkdir $N/bootstrap_correlation
fastspar_bootstrap --otu_table $input_file --number 1000 --prefix $N/bootstrap_counts/otu_data
done


# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
for item in `ls fastsparGTDB_in_141*/bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
parentdir2nd="$(dirname "$item")"
parentdir1st="$(dirname "$parentdir2nd")"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation $parentdir1st/bootstrap_correlation/cor_$N --covariance $parentdir1st/bootstrap_correlation/cov_$N
done

# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
for item in `ls fastsparGTDB_in_142*/bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
parentdir2nd="$(dirname "$item")"
parentdir1st="$(dirname "$parentdir2nd")"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation $parentdir1st/bootstrap_correlation/cor_$N --covariance $parentdir1st/bootstrap_correlation/cov_$N
done

# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
for item in `ls fastsparGTDB_in_143*/bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
parentdir2nd="$(dirname "$item")"
parentdir1st="$(dirname "$parentdir2nd")"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation $parentdir1st/bootstrap_correlation/cor_$N --covariance $parentdir1st/bootstrap_correlation/cov_$N
done

# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
for item in `ls fastsparGTDB_in_296*/bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
parentdir2nd="$(dirname "$item")"
parentdir1st="$(dirname "$parentdir2nd")"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation $parentdir1st/bootstrap_correlation/cor_$N --covariance $parentdir1st/bootstrap_correlation/cov_$N
done

# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
for item in `ls fastsparGTDB_in_297*/bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
parentdir2nd="$(dirname "$item")"
parentdir1st="$(dirname "$parentdir2nd")"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation $parentdir1st/bootstrap_correlation/cor_$N --covariance $parentdir1st/bootstrap_correlation/cov_$N
done

# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
for item in `ls fastsparGTDB_in_298*/bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
parentdir2nd="$(dirname "$item")"
parentdir1st="$(dirname "$parentdir2nd")"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation $parentdir1st/bootstrap_correlation/cor_$N --covariance $parentdir1st/bootstrap_correlation/cov_$N
done

# STEP 2: run the bootstrap (infer correlations for each bootstrap count)
for item in `ls fastsparGTDB_in_299*/bootstrap_counts/*`
do
filename=$(basename $item)     
N="${filename%.*}"
parentdir2nd="$(dirname "$item")"
parentdir1st="$(dirname "$parentdir2nd")"
fastspar -i 5 --threads 1 -y --otu_table $item --correlation $parentdir1st/bootstrap_correlation/cor_$N --covariance $parentdir1st/bootstrap_correlation/cov_$N
done


# STEP 3: run fastspar
for input_file in `ls fastsparGTDB_in_*.txt`
do
filename=$(basename $input_file)     
N="${filename%.*}"
fastspar -y --threads 1 --otu_table $input_file --correlation $N/median_correlation.tsv --covariance $N/median_covariance.tsv --iterations 5 
done


# STEP 4: fastspar pvalues (pvalues are calculated from the correlations)
for input_file in `ls fastsparGTDB_in_*.txt`
do
filename=$(basename $input_file)     
N="${filename%.*}"
fastspar_pvalues --threads 1 --otu_table $input_file --correlation $N/median_correlation.tsv --prefix $N/bootstrap_correlation/cor_ --permutations 1000 --outfile $N/pvalues.tsv
done


##########################################################################################################################################################
export your_dir=/shared/homes/s1/pig_microbiome/fastspar/fastspar_in_per_pig
qsub -V run_fastspar_per_pig.sh