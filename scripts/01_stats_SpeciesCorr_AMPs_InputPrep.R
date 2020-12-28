#This script requires the following packages:
install.packages("base", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("utils", repos = "http://cran.us.r-project.org")
install.packages("splitstackshape", repos = "http://cran.us.r-project.org")
install.packages("readr", repos = "http://cran.us.r-project.org")
install.packages("reshape2", repos = "http://cran.us.r-project.org")
install.packages("tidyr", repos = "http://cran.us.r-project.org")

rm(list = ls())

#upload all libraries
library(base)
library(readr) ## important
library(data.table)
library(dplyr) #impo
library(stringr)
library(utils)
library(reshape2) #important
library(tidyr)
library(splitstackshape) # important


#corr.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/macrel_stuff" # on local 
corr.dir = "/shared/homes/s1/pig_microbiome/fastspar/fastspar_in_per_pig" # on HPC

#macrel.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/macrel_stuff" # on local 
macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"


corr_data <- read.csv(file = file.path(corr.dir,"corr_pvalues_long_concat"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)

corr_data.edit <- corr_data 

sort(unique(corr_data.edit$speciesB))
# remove NA species
NROW(corr_data.edit)
corr_data.edit2 <- corr_data.edit[!grepl("NA__",corr_data.edit$speciesA),]
corr_data.edit2 <- corr_data.edit2[!grepl("NA__",corr_data.edit2$speciesB),]
NROW(corr_data.edit2)
corr_data.edit2$pig <- as.character(corr_data.edit2$pig)

corr_data.edit2$speciesA <- substr(corr_data.edit2$speciesA, 1, regexpr("\\__", corr_data.edit2$speciesA)-1)
corr_data.edit2$speciesB <- substr(corr_data.edit2$speciesB, 1, regexpr("\\__", corr_data.edit2$speciesB)-1)
corr_data.edit2$speciesA <- sub("^X", "", corr_data.edit2$speciesA)
corr_data.edit2$speciesB <- sub("^X", "", corr_data.edit2$speciesB)
corr_data.edit2$speciesA <- gsub("\\-", " ", corr_data.edit2$speciesA)
corr_data.edit2$speciesB <- gsub("\\.", " ", corr_data.edit2$speciesB)
head(corr_data.edit2)
tail(corr_data.edit2)

corr_data.edit2$species <- corr_data.edit2$speciesA

a <- sort(unique(corr_data.edit2$speciesA))
b <- sort(unique(corr_data.edit2$speciesB))
a <- as.data.frame(a)
b <- as.data.frame(b)
# if TRUE, same now match
NROW(a)==NROW(which(a==b))
###



###
merge_pig_contig_bin_species_clusters_out <- read_csv(file = file.path(macrel.dir,"merge_pig_contig_bin_species_clusters_out"),
                                                      col_types = cols(pig = col_character()))

amp_data <- as.data.frame(merge_pig_contig_bin_species_clusters_out)


### use the amp data to find which species are not present in the corr data, and remove them 
d <- sort(unique(amp_data$species))
d <- as.data.frame(d)
d$d <- gsub("\\-", " ", d$d)
head(d)
head(a)
colnames(d) <- "a"
# all the species from the corr dataframe that are not present in the AMP table: 
not_present <- as.list(anti_join(a, d, by = c("a")))

# remove them from the corr dataframe: 
corr_part1 <- subset(corr_data.edit2, !(speciesB %in% not_present$a))
corr_part2 <- subset(corr_part1, !(speciesA %in% not_present$a))
head(corr_data.edit2)
NROW(corr_data.edit2)
NROW(corr_part1)
NROW(corr_part2)
###

# continue with the amp data 
amp_data$species <- gsub("\\-", " ", amp_data$species)



test <- amp_data %>% 
  dplyr::select(pig, species, cluster_98) %>%
  group_by(pig,cluster_98) %>% 
  tally() %>%
  dplyr::arrange(desc(n))

test2 <- amp_data %>% 
  dplyr::select(pig, species, cluster_98) %>%
  group_by(cluster_98) %>% 
  tally() %>%
  dplyr::arrange(desc(n))

#############################################
#############################################
## IMPORTANT TO RUN ON LOCAL
#############################################
#############################################
# amp_data <- amp_data %>%
#   dplyr::filter(cluster_98=="Cluster 22071"|
#                   cluster_98=="Cluster 15993"|
#                   cluster_98=="Cluster 11017"|
#                   cluster_98=="Cluster 25291"|
#                   cluster_98=="Cluster 19756")
#############################################
#############################################

amp_data$pig_species <- paste0(amp_data$pig,".",amp_data$species)
amp_data_sub <- amp_data %>% 
  dplyr::select(pig_species, cluster_98) 


si <- unique(amp_data_sub)
NROW(si)
# remove species that contain NA (these contigs don't associate to any species so it's worthless taking them on to look at species corr)
si <- si[!grepl(".NA",si$pig_species),]
NROW(si)


amp_data.binary <- reshape2::dcast(si, formula = pig_species ~ cluster_98, fun.aggregate = length)

# split column 
amp_data.binary2 <- cSplit(amp_data.binary, "pig_species", ".")
# move last two cols ahead
amp_data.binary3 <- amp_data.binary2 %>%
  dplyr::select(pig_species_1, pig_species_2, everything()) 
names(amp_data.binary3)[names(amp_data.binary3) == 'pig_species_1'] <- 'pig'
names(amp_data.binary3)[names(amp_data.binary3) == 'pig_species_2'] <- 'species'
amp_data.binary3$pig <- as.character(amp_data.binary3$pig)
##


# save the results 
fwrite(x=amp_data.binary3, file = file.path(macrel.dir,"amp_data.binary_clean"))

# save the results 
fwrite(x=corr_data.edit2, file = file.path(macrel.dir,"corr_data_clean"))

#View(amp_data.binary3)
#View(corr_data.edit2)
