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


# open the (mean) correlations files 
pea <- read_delim(paste0(out_dir,"species_corr_pearson_mean.tsv"), 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
sp <- read_delim(paste0(out_dir,"species_corr_spearman_mean.tsv"), 
                 "\t", escape_double = FALSE, trim_ws = TRUE)


pea$method <- paste0("Pearson")
sp$method <- paste0("Spearman")

pea_sp <- rbind(pea,sp)

#########


# VISUALIZE (mean) correlations via heatmaps: 


return_matrix_for_heatmap <- function (df) {
  
  df1 <- cSplit(df, "pairs",".")
  
  df2 <- df1 %>% 
    dplyr::select(pairs_1, pairs_2, mean_corr) %>% 
    group_by(pairs_1) %>% 
    arrange()
  
  df3 <- df2 %>% 
    pivot_wider(names_from=pairs_1,values_from=mean_corr, values_fill = 0) %>% 
    group_by(pairs_2) %>% 
    arrange()
  
  df3 <- as.data.frame(df3)
  
  rownames(df3) <- df3$pairs_2
  df3$pairs_2 <- NULL
  data <- as.matrix(df3)
  
  return(data)
}




# HEATMAP of pearson correlations
d1 <- return_matrix_for_heatmap(pea)
dim(d1)
# subset it
d1 <- d1[300:400,300:400]
pdf(paste0(out_dir, "heatmap_pearson.pdf"))
pheatmap(d1, main = "no_clu", cluster_cols = F, cluster_rows = F, fontsize = 2)
pheatmap(d1, main = "clustered", cluster_cols = T, cluster_rows = T, fontsize = 2)
dev.off()


# HEATMAP of spearman correlations
d2 <- return_matrix_for_heatmap(sp)
dim(d2)
# subset it
d2 <- d2[300:400,300:400]
pdf(paste0(out_dir, "heatmap_spearman.pdf"))
pheatmap(d2, main = "no_clu", cluster_cols = F, cluster_rows = F, fontsize = 2)
pheatmap(d2, main = "clustered", cluster_cols = T, cluster_rows = T, fontsize = 2)
dev.off()



#################################################################################
#################################################################################

# interactive heatmap
d1 <- return_matrix_for_heatmap(pea)
dim(d1)

# subset it
d1_mini <- d1[1:100,1:100]

heatmaply(d1_mini, cexRow = 0.1, cexCol = 0.1,
          show_dendrogram = FALSE, 
          file = paste0(out_dir, "heatmap_interactive.html"))

# # subset it
# d2_mini <- d2[300:350,300:350]
# heatmaply(d2_mini, cexRow = 0.5, cexCol = 0.5,
#           show_dendrogram = FALSE)


#################################################################################
#################################################################################
