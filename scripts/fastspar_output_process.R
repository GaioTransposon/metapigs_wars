library(readr)
library(pheatmap)
library(reshape2)


out_dir = "/Users/12705859/Desktop/metapigs_wars/out/"  # local



median_correlation_editUnderscore <- read_table2("~/Desktop/metapigs_dry/network_analysis/fastspar/bigpigs_outtie/median_correlation_editUnderscore.tsv")


median_correlation_editUnderscore <- as.data.frame(median_correlation_editUnderscore)
rownames(median_correlation_editUnderscore) <- median_correlation_editUnderscore[,1]
median_correlation_editUnderscore[,1] <- NULL
colnames(median_correlation_editUnderscore)



# melt to look at distribution of correlations:

dat <- melt(as.matrix(median_correlation_editUnderscore))
hist(dat$value)
closeAllConnections()
dat %>% 
  dplyr::filter(Var1=="557_1") %>%
  dplyr::filter(Var2=="354_0")


# subset it
test <- median_correlation_editUnderscore[350:450,350:450]

# HEATMAP
pdf(paste0(out_dir, "heatmap_fastspar.pdf"))
pheatmap(test, main = "no_clu", cluster_cols = F, cluster_rows = F, fontsize = 2)
pheatmap(test, main = "clustered", cluster_cols = T, cluster_rows = T, fontsize = 5)
dev.off()



# # HEATMAP all
# pdf(paste0(out_dir, "heatmap_fastspar_large.pdf"))
# pheatmap(median_correlation_editUnderscore, main = "no_clu", cluster_cols = F, cluster_rows = F, fontsize = 2)
# pheatmap(median_correlation_editUnderscore, main = "clustered", cluster_cols = T, cluster_rows = T, fontsize = 2)
# dev.off()

