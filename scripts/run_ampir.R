
# download :
#install.packages("Matrix", repos = "http://cran.us.r-project.org", dependencies=TRUE)
#install.packages("ampir", repos = "http://cran.us.r-project.org", dependencies=TRUE)

#upload all libraries
library("base")
library("utils")
library("ampir")

input_file_dir = "/shared/homes/12705859/prodigal_work/"
out_dir = "/shared/homes/12705859/ampir/"

# test was done on mini input file: mini_noAsteriks.faa
# produced from `head all_concatenated_AA_noAsteriks.faa -n 20000 > mini_noAsteriks.faa`

# read in input:
my_protein_df <- read_faa(file = paste0(input_file_dir,"mini_noAsteriks.faa"))    
# when testing on mini dataset:
#my_protein_df <- read_faa(file = paste0(input_file_dir,"all_concatenated_AA_noAsteriks.faa")) 

# Calculate the probability that each protein is an antimicrobial peptide with predict_amps(). 
my_prediction <- predict_amps(my_protein_df, model = "precursor")

# Predicted proteins with a specified predicted probability value could then be extracted and written to a FASTA file:
my_predicted_amps <- my_protein_df[which(my_prediction$prob_AMP >= 0.8),]

# Write the data.frame with sequence names in the first column and protein sequences in the second column to a FASTA formatted file with df_to_faa()
df_to_faa(my_predicted_amps, paste0(out_dir,"my_predicted_amps.fasta"))


####
# Input file was obtained from prodigal_on_assemblies.sh, then * removed with: 
# sed 's|[*]||g' all_concatenated_assemblies_AA.faa > all_concatenated_AA_noAsteriks.faa
####