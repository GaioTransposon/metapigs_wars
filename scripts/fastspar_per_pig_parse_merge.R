# runs from the HPC 
# language: R 

#This script requires the following packages:
install.packages("base", repos = "http://cran.us.r-project.org")
install.packages("readr", repos = "http://cran.us.r-project.org")
#install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE) # usually already installed
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("utils", repos = "http://cran.us.r-project.org")
install.packages("reshape", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(readr)
library(data.table)
library(dplyr)
library(stringr)
library(utils)
library(reshape)


#corr.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs" # on local 
corr.dir = "/shared/homes/s1/pig_microbiome/fastspar/fastspar_in_per_pig" # on HPC


dirs <- list.dirs(corr.dir, recursive = FALSE)
dirs <- dirs[ grepl("fastsparGTDB_in_", dirs) ]

for (corr_dir in dirs) {
  
  raw_median_corr_files = list.files(corr_dir, pattern="median_correlation.tsv")
  
  df <- read.csv(file.path(corr_dir,raw_median_corr_files), sep = "\t", 
                     row.names = NULL, header = TRUE, stringsAsFactors = FALSE)

  df <- as.data.frame(df)
  rownames(df) <- df[,1]
  df[,1] <- NULL
  
  dat_corr <- melt(as.matrix(df))
  colnames(dat_corr) <- c("speciesA","speciesB","median_corr")
  ###
  
  raw_pvalues_files = list.files(corr_dir, pattern="pvalues.tsv")
  
  df1 <- read.csv(file.path(corr_dir,raw_pvalues_files), sep = "\t", 
                 row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
  
  df1 <- as.data.frame(df1)
  rownames(df1) <- df1[,1]
  df1[,1] <- NULL
  
  dat_pvalues <- melt(as.matrix(df1))
  colnames(dat_pvalues) <- c("speciesA","speciesB","pvalue")
  
  ###
  
  dat <- full_join(dat_corr,dat_pvalues)
  dat$pig <- paste0(as.character(str_remove(basename(corr_dir),"fastsparGTDB_in_")))
    
  # save data 
  fwrite(x = dat,
         file = file.path(corr_dir,"corr_pvalues_long"))

}


# merge all 

df0 = data.frame(
  speciesA = character(),
  speciesB = character(),
  median_corr = numeric(),
  pvalues = numeric(),
  pig = character()
)

for (corr_dir in dirs) {
  
  long_files = list.files(corr_dir, pattern="corr_pvalues_long")
  
  df <- read.csv(file.path(corr_dir,long_files), sep = ",", 
                 row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
  
  df1 <- as.data.frame(df)
  
  df0 <- rbind(df0, df1)
}

head(df0)
tail(df0)

# save data 
fwrite(x = df0,
       file = file.path(corr.dir,"corr_pvalues_long_concat"))

# df1 <- df0 %>% 
#   dplyr::filter(pig=="fastsparGTDB_in_test1")
# View(df1)
# hist(df1$median_corr)
# hist(df1$pvalue)
# 
# 
# library(ggpubr)
# library(rstatix)
# library(emmeans)
# df1$pig <- as.factor(df1$pig)
# res.aov <- df1 %>% anova_test(median_corr ~ pig)
# res.aov


     