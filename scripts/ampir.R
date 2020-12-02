#!/usr/bin/Rscript
# version this script was developed in: "R version 3.6.0 (2019-04-26)"
# platform this script was developed in: "x86_64-apple-darwin15.6.0"


# STEP 1: create new R env on HPC (doing on screen: 1292879.newRenv) --> conda create -n R361_env -c conda-forge r-base=3.6.1
# STEP 2: conda activate R361_env, and install devtools: conda install -c conda-forge r-devtools
# STEP 3: enter R, then install.packages("cli"); then install.packages("rlang", type = "source"), then exit and reenter R
# STEP 4: devtools::install_github("Legana/ampir") (say yes to update all packages)
# asked me to remove a directory, done, re-do step 4. 
# again a problem: namespace ‘glue’ 1.3.1 is already loaded, but >= 1.3.2 is required
# so exit, re-enter, then: install.packages("glue", type = "source"), re-do step 4. didn't work
#STEP 4: install.packages("ampir", type = "source") (worked)
# STEP 4: run the script below

screen -S run_ampir
conda activate R361_env
install.packages("ampir", type = "source")



#upload all libraries
library(base)
library(utils)
library(ampir)

# input dir
input_file_dir = "/shared/homes/12705859/prodigal_work/"
out_dir = "/shared/homes/12705859/ampir/"


# test file: 
# head -n 100000 all_concatenated_AA.faa > mini
# sed 's|[*]||g' mini > mini2

# pre-processing on input file: 
# sed 's|[*]||g' all_concatenated_AA.faa > all_concatenated_AA_noAsteriks.faa

# read in 
my_protein_df <- read_faa(file = paste0(input_file_dir,"all_concatenated_AA_noAsteriks.faa"))

# Calculate the probability that each protein is an antimicrobial peptide with predict_amps(). 
# Since these proteins are all full length precursors rather than mature peptides we use ampir’s built-in precursor model.
# Note that amino acid sequences that are shorter than 10 amino acids long and/or contain anything other than the standard 20 amino acids are not evaluated and will contain an NA as their prob_AMP value.
#  The default model, “precursor” is best suited for full length proteins and the “mature” model is best suited for small mature proteins (<60 amino acids)
# so precursos is fine for the AA predicted from my bins 
my_prediction <- predict_amps(my_protein_df, model = "precursor")

# Predicted proteins with a specified predicted probability value could then be extracted and written to a FASTA file:
my_predicted_amps <- my_protein_df[which(my_prediction$prob_AMP >= 0.8),]

# Write the data.frame with sequence names in the first column and protein sequences in the second column to a FASTA formatted file with df_to_faa()
df_to_faa(my_predicted_amps, paste0(out_dir,"my_predicted_amps.fasta"))




conda activate R361_env
cd ampir

