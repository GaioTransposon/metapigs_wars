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
# dictionary

# OUTPUTS:

# lots!!


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



# splitting into multiple dataframes (by pigID)
multiple_DFs <- split( df3 , f = df3$pig )


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
#################################################################################








#########

# dplyr summarise to get the mean, median, and sd 
# store this info 

final_df$pairs <- paste0(final_df$Var1,".",final_df$Var2)
  
# reduce to keep: mean, Var1, Var2 (the two species) 
df <- final_df %>% 
  dplyr::select(pairs,corr_value) %>% 
  group_by(pairs) %>% 
  dplyr::summarise(mean_corr=mean(corr_value),
                   median_corr=median(corr_value),
                   sd_corr=sd(corr_value),
                   subjects=n())

head(df)
View(df)
hist(df$mean_corr)

# store this info (fwrite) <- TO DO !

test <- df

# sd NA are the correlatrions found in only one subject, so we can discard those. 
test <- test[complete.cases(test), ]

# pos correlations
test_pos <- test %>%
  dplyr::filter(sd_corr < 0.1) %>% 
  dplyr::filter(mean_corr > 0.5) 

# neg correlations
test_neg <- test %>%
  dplyr::filter(sd_corr < 0.1) %>% 
  dplyr::filter(mean_corr < 0.5) 




# filter by bug of interest
# let's say you want to know what probiotics can be effective to suppress the growth of E. coli
test_Esc <- test_neg[test_neg$pairs %like% "Escher", ]
View(test_Esc)


#########

# show top (most often occurring) negatively correlating bugs
test_neg %>% 
  arrange(desc(n)) 




head(final_df)
list_pigs <- final_df %>% 
  dplyr::filter(pairs=="Cloacibacillus porcorum__11272.Prevotella sp000436035__12268") %>% 
  dplyr::filter(corr_value < 0.5) %>%
  dplyr::select(pig)
list_pigs <- as.list(unique(list_pigs$pig))

# let's visualize the time trend of two pos (or neg) correlating bugs from the raw data: 
W <- df3
W <- subset(W, (pig %in% list_pigs))
head(W)

pdf("~/Desktop/cors.pdf")
W %>% 
  dplyr::filter(gOTU=="Prevotella sp000436035__12268" | 
                  gOTU=="Cloacibacillus porcorum__11272") %>% 
  ggplot(., aes(x=date,y=log(norm_value), group=gOTU, color=gOTU)) +
  geom_line() + geom_point(size=0.8)+
  theme_bw() +
  theme(axis.text.y = element_text(size=1)) +
  facet_wrap(~pig, scales="free_y")
dev.off()




#########

# VISUALIZE via heatmaps: 


#df1 <- cSplit(df, "pairs",".")

df2 <- df1 %>% 
  dplyr::select(pairs_1, pairs_2, mean_corr) %>% 
  group_by(pairs_1) %>% 
  arrange()
NROW(df2)

df3 <- df2 %>% 
  pivot_wider(names_from=pairs_1,values_from=mean_corr, values_fill = 0) %>% 
  group_by(pairs_2) %>% 
  arrange()

df3 <- as.data.frame(df3)
head(df3)

rownames(df3) <- df3$pairs_2
df3$pairs_2 <- NULL
data <- as.matrix(df3)
dim(data)

# subsets
data <- data[300:400,300:400]
#data <- data[1:10,1:10]


pdf("~/Desktop/friends_and_enemies.pdf")
pheatmap(data, cluster_cols = F, cluster_rows = F, fontsize = 2)
pheatmap(data, cluster_cols = T, cluster_rows = T, fontsize = 2)
dev.off()


