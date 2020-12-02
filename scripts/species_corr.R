library(corrplot)
library(RColorBrewer)
library(readr)
library(splitstackshape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(data.table)
library(pheatmap)
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

colnames(df0)[colnames(df0) == 'node'] <- 'gOTU'

df0$gOTU <- as.character(df0$gOTU)

df0$gOTU <- paste0(df0$species,"__",df0$gOTU)

# cohort selection (all piggies, no mothers, or pos/neg controls)
df0 <- df0 %>% dplyr::filter(cohort=="Control" | 
                               cohort=="D-Scour" | 
                               cohort=="ColiGuard" | 
                               cohort=="Neomycin" | 
                               cohort=="NeoD" | 
                               cohort=="NeoC" )

######################################################################

# CREATE COUNTS TABLE , NORMALIZE by LIB SIZE
df1 <- df0
head(df1)
NROW(df1)

### minitest: 
# filter a sample (pig,date)
test <- df1 %>% dplyr::filter(pig=="14159") %>% dplyr::filter(date=="t2")
sum(test$value)
NROW(test)
# for each sample (pig,date), sum up together the counts that fall within one species (same species assigned to distinct bins)
test2 <- test %>%
  group_by(pig,gOTU,date) %>%
  dplyr::summarize(sum_value = sum(value)) 
sum(test2$sum_value)
NROW(test2)
# normalize by library size 
test3 <- test2 %>% 
  group_by(pig,date) %>% 
  dplyr::mutate(norm_value = sum_value/sum(sum_value)) %>% 
  dplyr::select(-sum_value)
sum(test3$norm_value)
NROW(test3)
###


# PROCEED to all: 

# for each sample (pig,date), sum up the counts that fall within one species (same species assigned to distinct bins)
df2 <- df1 %>%
  group_by(pig,gOTU,date) %>%
  dplyr::summarize(sum_value = sum(value)) 
head(df2)
sum(df2$sum_value)

# normalize by library size 
df3 <- df2 %>% 
  group_by(pig,date) %>% 
  dplyr::mutate(norm_value = sum_value/sum(sum_value)) %>% 
  dplyr::select(-sum_value)
head(df3)

# if your total sum is equal to the total number of samples, 
# it means that the sum within each sample (pig,date) is 1, and that's correct  
NROW(unique(paste0(df3$pig,df3$date)))==sum(df3$norm_value)


df3 <- as.data.frame(df3)
# save normalized data
fwrite(df3, file=paste0(out_dir,"no_reps_all_norm.tsv"), sep = "\t")


# splitting into multiple dataframes (by pigID)
multiple_DFs <- split( df3 , f = df3$pig, drop = TRUE)

#################################################################################

# Pearson correlation

# construct an empty dataframe to build on 
final_df_pea <- data.frame(
  Var1 = character(),
  Var2 = character(),
  corr_value = character(),
  pig = numeric(),
  stringsAsFactors = FALSE
)

for (single_DF in multiple_DFs) {
  
  # correlation is computed for subjects with more than 1 time point (obviously)
  if ( NROW(unique(single_DF$date)) > 1 ) {
    
    single_DF <- as.data.frame(single_DF)
    colnames(single_DF) <- c("pig","gOTU","date","norm_value")
    
    pigID <- unique(single_DF$pig)
    
    # start 
    test <- single_DF %>% 
      dplyr::select(gOTU,date,norm_value) 
    
    test <- test %>% 
      pivot_wider(names_from=gOTU,values_from=norm_value)
    
    test <- as.data.frame(test)
    
    rownames(test) <- test$date
    test$date <- NULL
    
    M <-cor(test, method = "pearson")
    
    #####
    # if you want you can plot the single corr plot: (must ".pdf" it)
    # corrplot(M, type="upper", order="hclust",
    #          col=brewer.pal(n=8, name="RdYlBu"))
    #####
    
    # convert to df: 
    
    df <- na.omit(melt(M))  # reshaping
    df <- df[order(df$Var1), ]   # ordering
    
    # add column with pig ID
    # end of loop
    df$pig <- paste0(as.character(pigID))
    
    colnames(df) <- c("Var1", "Var2", "corr_value","pig") # setting colnames
    
    # rbind with all other pigs
    final_df_pea <- rbind(
      final_df_pea,
      df
    )
    
  } else {
    
    print("Subject with only one time point sample excluded")
  }
    
}


head(final_df_pea)
tail(final_df_pea)
NROW(unique(df0$pig))
hist(final_df_pea$corr_value)

# save file
fwrite(final_df_pea, file=paste0(out_dir,"species_corr_pearson.tsv"), sep = "\t")


#################################################################################

# Spearman correlation


# construct an empty dataframe to build on 
final_df_sp <- data.frame(
  Var1 = character(),
  Var2 = character(),
  corr_value = character(),
  pig = numeric(),
  stringsAsFactors = FALSE
)

for (single_DF in multiple_DFs) {
  
  # correlation is computed for subjects with more than 1 time point (obviously)
  if ( NROW(unique(single_DF$date)) > 1 ) {
    
    single_DF <- as.data.frame(single_DF)
    colnames(single_DF) <- c("pig","gOTU","date","norm_value")
    
    pigID <- unique(single_DF$pig)
    
    # start 
    test <- single_DF %>% 
      dplyr::select(gOTU,date,norm_value) 
    
    test <- test %>% 
      pivot_wider(names_from=gOTU,values_from=norm_value)
    
    test <- as.data.frame(test)
    
    rownames(test) <- test$date
    test$date <- NULL
    
    M <-cor(test, method = "spearman")
    
    #####
    # if you want you can plot the single corr plot: (must ".pdf" it)
    # corrplot(M, type="upper", order="hclust",
    #          col=brewer.pal(n=8, name="RdYlBu"))
    #####
    
    # convert to df: 
    
    df <- na.omit(melt(M))  # reshaping
    df <- df[order(df$Var1), ]   # ordering
    
    # add column with pig ID
    # end of loop
    df$pig <- paste0(as.character(pigID))
    
    colnames(df) <- c("Var1", "Var2", "corr_value","pig") # setting colnames
    
    # rbind with all other pigs
    final_df_sp <- rbind(
      final_df_sp,
      df
    )
    
  } else {
    
    print("Subject with only one time point sample excluded")
  }
  
}


# save file
fwrite(final_df_sp, file=paste0(out_dir,"species_corr_spearman.tsv"), sep = "\t")


#################################################################################

# transform to mean, median, sd of correlations and save


#####

final_df_pea$pairs <- paste0(final_df_pea$Var1,".",final_df_pea$Var2)

# reduce to keep: mean, Var1, Var2 
pea <- final_df_pea %>% 
  dplyr::select(pairs,corr_value) %>% 
  group_by(pairs) %>% 
  dplyr::summarise(mean_corr=mean(corr_value),
                   median_corr=median(corr_value),
                   sd_corr=sd(corr_value),
                   n_subjects=n())

pdf(paste0(out_dir,"histogram_pea_corr.pdf"))
hist(pea$mean_corr)
dev.off()

# save file
fwrite(pea, file=paste0(out_dir,"species_corr_pearson_mean.tsv"), sep = "\t")

#####


final_df_sp$pairs <- paste0(final_df_sp$Var1,".",final_df_sp$Var2)

# reduce to keep: mean, Var1, Var2 
sp <- final_df_sp %>% 
  dplyr::select(pairs,corr_value) %>% 
  group_by(pairs) %>% 
  dplyr::summarise(mean_corr=mean(corr_value),
                   median_corr=median(corr_value),
                   sd_corr=sd(corr_value),
                   n_subjects=n())

pdf(paste0(out_dir,"histogram_pea_corr.pdf"))
hist(sp$mean_corr)
dev.off()

# save file
fwrite(sp, file=paste0(out_dir,"species_corr_spearman_mean.tsv"), sep = "\t")


#################################################################################
#################################################################################

