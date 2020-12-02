
# dRep_output_analysis.R
# analysis of dRep output 

library(readxl)
library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(splitstackshape)
library(pheatmap)
library(ggpubr)
library(robCompositions)
library(ggbiplot)

source_dir = "/Users/12705859/metapigs_dry/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_dry/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_dry/dRep/"  # local


# upload dRep output 
Cdb <- read_csv(paste0(middle_dir,"Cdb.csv"))

# upload cohorts info
cohorts <- read_xlsx(paste0(source_dir,"cohorts.xlsx"))

C1 <- separate(data = Cdb, col = genome, into = c("pig", "bin"), sep = "_")
C1 <- C1[,c("pig","bin","primary_cluster","secondary_cluster")]
C1$primary_cluster <- as.character(C1$primary_cluster)
head(C1)

# upload bins with counts (from output of 7.R)
no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)

# remove .fa extension to match bins in checkm df 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

no_reps_all$primary_cluster <- paste0(no_reps_all$secondary_cluster)
no_reps_all <- cSplit(no_reps_all,"primary_cluster","_")
no_reps_all$primary_cluster_2 <- NULL
colnames(no_reps_all)[colnames(no_reps_all)=="primary_cluster_1"] <- "primary_cluster"

######################################################################

# load gtdbtk assignments of the bins

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(middle_dir,"gtdb_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


######################################################################


# create text file to contain dRep text output

sink(file = paste0(out_dir,"dRep_numbers.txt"), 
     append = FALSE, type = c("output"))
sink()


###########################


sink(file = paste0(out_dir,"dRep_numbers.txt"), 
     append = TRUE, type = c("output"))
paste0("dRep-clustered bins: ", 
       round(NROW(C1)/NROW(gtdbtk_bins)*100,2),
       "%",
       " (n=",NROW(C1),")" )
paste0("of which primary clusters: ", length(unique(C1$primary_cluster)) )
paste0("of which secondary clusters ", length(unique(C1$secondary_cluster)) )
sink()



########################################################################################################


# Extent of agreement between dRep and GTDBTK classification: 


df <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))


NROW(df)
z <- df %>% 
  dplyr::select(pig,bin,secondary_cluster,domain,phylum,class,order,family,genus,species,node) %>%
  distinct() %>% 
  dplyr::filter(!secondary_cluster=="no_cluster")
NROW(z)

z$digit <- 1

zz <- z %>% group_by(species,secondary_cluster) %>%
  dplyr::summarize(sum_digit=sum(digit)) %>% 
  group_by(secondary_cluster) %>% 
  dplyr::mutate(freq=sum_digit/sum(sum_digit))

View(zz)
