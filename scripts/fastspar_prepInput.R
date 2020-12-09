library(RColorBrewer)
library(readr)
library(splitstackshape)
library(tidyr)
library(dplyr)
library(matrixStats)
library(data.table)
library(readxl)
library(ggpubr)
library(forcats)
library(data.table)
library(reshape2)

source_dir = "/Users/12705859/metapigs_wars/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_wars/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_wars/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_wars/out/"  # local

######################################################################

# input files: 
# gtdbtk_bins_completeTaxa
# no_reps_all.csv (BINS COUNTS)

# OUTPUTS:

# species_corr_pearson.tsv
# species_corr_spearman.tsv
# species_corr_pearson_mean.tsv
# species_corr_spearman_mean.tsv


######################################################################

# counts data 

no_reps_all <- read.csv(paste0(middle_dir,"no_reps_all.csv"), 
                        na.strings=c("","NA"),
                        check.names = FALSE,
                        header = TRUE)


# remove .fa extension 
no_reps_all$bin <- gsub(".fa","", no_reps_all$bin)
head(no_reps_all)
NROW(no_reps_all)

######################################################################

# load gtdbtk assignments of the bins
gtdbtk_bins <- read_csv(paste0(middle_dir,"gtdb_bins_completeTaxa"),
                        col_types = cols(node = col_character(),
                                         pig = col_character()))


head(gtdbtk_bins)

######################################################################

# merge bins info to gtdbtk assignment info :  


df0 <- merge(no_reps_all, gtdbtk_bins, by=c("pig","bin"))

# rename node as gOTU and place "gOTU_" in front of node number: a separate genomic OTU identifier for each different genome

colnames(df0)[colnames(df0) == 'node'] <- '#OTU ID'

df0$`#OTU ID` <- as.character(df0$`#OTU ID`)

df0$`#OTU ID` <- paste0(df0$species,"__",df0$`#OTU ID`)

# cohort selection (all piggies, no mothers, or pos/neg controls)
df0 <- df0 %>% dplyr::filter(cohort=="Control" | 
                               cohort=="D-Scour" | 
                               cohort=="ColiGuard" | 
                               cohort=="Neomycin" | 
                               cohort=="NeoD" | 
                               cohort=="NeoC" )

######################################################################


df1 <- df0

# the only transform we need to do is to sum up the counts that fall within one species (same species assigned to distinct bins)
# this should not happen but it does happen. 
# Either because some contigs belonging to the same organism don t get binned together or because of taxa mis-classification. 
df2 <- df1 %>%
  group_by(pig,`#OTU ID`,date) %>%
  dplyr::summarize(sum_value = sum(value)) 
head(df2)



z <- df2 %>% 
  dplyr::select(sum_value,`#OTU ID`,pig,date)


z$sample <- paste0(z$date,"_",z$pig)

z <- as.data.frame(z)


z <- z %>% 
  dplyr::select(sum_value,`#OTU ID`,sample) %>% 
  dplyr::mutate(sum_value=round(sum_value)) %>% 
  pivot_wider(names_from=sample, values_from=sum_value, values_fill=0)



z$`#OTU ID`
colnames(z)

write.table(x = z, file = "~/Desktop/metapigs_wars/fastsparGTDB_in.txt", row.names = FALSE, quote = FALSE, sep = '\t')
colnames(fastspar_Mothers_in)

