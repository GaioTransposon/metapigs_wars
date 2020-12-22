# /shared/homes/12705859/contig_abundances/no_reps_contigs_PigBinContig.csv 1762877983
# pig,bin,contigName
# 29667,bins.13.fa,k141_100
# 29778,bins.130.fa,k141_100
# head no_reps_contigs_PigBinContig.csv -n 1000 > test1
# 
# /shared/homes/12705859/cdhit_work/cd_hit_onMacrel/contigs_all_clusters 5964739
# pig,contig,ORF_number,cluster_100,cluster_90,cluster_95,cluster_98
# 14171,k141_189950,1,Cluster 0,Cluster 0,Cluster 0,Cluster 0
# 14200,k141_40977,1,Cluster 0,Cluster 0,Cluster 0,Cluster 0
# cp contigs_all_clusters test2
# 
# 
# /shared/homes/s1/pig_microbiome/kraken_on_concat_bins/all_concatenated_220_essential.kraken 3250400
# C	14159_bins.100.fa	CAG-914 sp000437895 (taxid 15578)	1230079
# C	14159_bins.101.fa	PeH17 sp001940845 (taxid 25959)	1543424
# cp all_concatenated_220_essential.kraken test3


# runs from the HPC 
# language: R 

#This script requires the following packages:
install.packages("base", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("utils", repos = "http://cran.us.r-project.org")
install.packages("splitstackshape", repos = "http://cran.us.r-project.org")
install.packages("readr", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(readr)
library(data.table)
library(dplyr)
library(stringr)
library(utils)
library(splitstackshape)


macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"

fileA <- read_csv("/shared/homes/12705859/contig_abundances/no_reps_contigs_PigBinContig.csv", col_types = cols(pig = col_character()))

fileB <- read.csv(file.path(macrel.dir,"contigs_all_clusters"), sep = ",", 
                  row.names = NULL, header = TRUE, stringsAsFactors = FALSE)

fileC <- read.csv("/shared/homes/s1/pig_microbiome/kraken_on_concat_bins/all_concatenated_220_essential.kraken", sep = "\t", 
                  row.names = NULL, header = FALSE, stringsAsFactors = FALSE)

## run on local:
# fileA <- read_csv("Desktop/bins_clustering_parsing_DFs/parsing_metapigs_wars/test1", col_types = cols(pig = col_character()))
# 
# fileB <- read.csv("Desktop/bins_clustering_parsing_DFs/parsing_metapigs_wars/test2", sep = ",", 
#                   row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
# 
# fileC <- read.csv("Desktop/bins_clustering_parsing_DFs/parsing_metapigs_wars/test3", sep = "\t", 
#                   row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
# 

colnames(fileA)[colnames(fileA)=="contigName"] <- "contig"


fileC <- cSplit(fileC, "V2", "_")
fileC <- cSplit(fileC, "V3", "(")

fileC <- fileC %>%
  dplyr::select(V2_1, V2_2, V3_1)

colnames(fileC) <- c("pig","bin","species")
head(fileC)

head(fileA)
head(fileB)

head(fileC)

NROW(fileA)
NROW(fileB)
fileAB <- right_join(fileA,fileB)

NROW(fileAB)
NROW(fileC)
fileABC <- right_join(fileC,fileAB)
NROW(fileABC)
head(fileABC)

fileABC <- fileABC %>% 
  dplyr::select(species, everything()) %>%
  dplyr::select(bin, everything()) %>%
  dplyr::select(contig, everything()) %>%
  dplyr::select(pig, everything())


# save 
fwrite(x = fileABC,
       file = file.path(macrel.dir,"merge_pig_contig_bin_species_clusters_out"))

