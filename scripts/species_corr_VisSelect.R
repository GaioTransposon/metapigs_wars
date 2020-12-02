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
library(heatmaply)



source_dir = "/Users/12705859/metapigs_wars/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_wars/middle_dir/" # git 
out_dir_git = "/Users/12705859/metapigs_wars/out/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_wars/out/"  # local

######################################################################

# input files: 
# species_corr_pearson_mean.tsv
# species_corr_spearman_mean.tsv

# OUTPUTS:

# .....


######################################################################

# open the (mean) correlations files 
no_reps_all_norm <- no_reps_all_norm <- read_delim(paste0(out_dir, "no_reps_all_norm.tsv"), 
                                                   "\t", escape_double = FALSE, col_types = cols(pig = col_character()), 
                                                   trim_ws = TRUE)


# open the correlations files 
pea <- read_delim(paste0(out_dir,"species_corr_pearson.tsv"), 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
sp <- read_delim(paste0(out_dir,"species_corr_spearman.tsv"), 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
pea <- as.data.frame(pea)
sp <- as.data.frame(sp)

# open the mean correlations files 
pea_mean <- read_delim(paste0(out_dir,"species_corr_pearson_mean.tsv"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
sp_mean <- read_delim(paste0(out_dir,"species_corr_spearman_mean.tsv"), 
                 "\t", escape_double = FALSE, trim_ws = TRUE)



#########

# let's continue with pearson only: 

# sd NA are the correlations found in only one subject, so we can discard those. 
NROW(pea_mean)
pea_mean <- pea_mean[complete.cases(pea_mean), ]
NROW(pea_mean)
#########


# pos correlations
pea_pos <- pea_mean %>%
  dplyr::filter(sd_corr < 0.1) %>% 
  dplyr::filter(mean_corr > 0.5) %>% 
  arrange(desc(n_subjects)) 

# neg correlations
pea_neg <- pea_mean %>%
  dplyr::filter(sd_corr < 0.1) %>% 
  dplyr::filter(mean_corr < -0.3) %>% 
  arrange(desc(n_subjects)) %>% 
  top_n(20) 
pairs_list <- as.list(unique(pea_neg$pairs))
pea_neg <- cSplit(pea_neg, "pairs",".")
colnames(pea_neg) <- c("mean_corr","median_corr","sd_corr","n_subjects", "Var1","Var2")
head(pea_neg)

# based on these negative correlations we want to find
# which of the piglets showed these negative correlations: 
NROW(pea_neg)
NROW(pea)
merged <- right_join(pea,pea_neg)
NROW(merged)
head(merged)

merged2 <- merged %>% 
  dplyr::filter(corr_value < -0.4)
pig_list <- as.list(unique(merged2$pig))
NROW(pig_list)


# and show the abundance across time from the normalized abundance dataframe


# loop version 2
for (i in 1:length(pairs_list)) {
  n <- print(pairs_list[[i]])
}

n


NROW(new)
head(no_reps_all_norm)
sub <- no_reps_all_norm %>%
  dplyr::filter(gOTU == "Cloacibacillus porcorum__11272" |
                  gOTU == "Prevotella sp000436035__12268")
sub <- subset(sub, (pig %in% mylist))
pdf(paste0(out_dir,"neg_corr_species.pdf"))
sub %>% 
  ggplot(., aes(x=date,y=log(norm_value),group=gOTU,color=gOTU))+
  geom_line()+ 
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(size=0.2),
        axis.text.y = element_text(size=0.2),
        legend.position = "top")+
  facet_wrap(~pig, scales = "free_y")
dev.off()

