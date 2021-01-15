#This script requires the following packages:
install.packages("base", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("utils", repos = "http://cran.us.r-project.org")
install.packages("readr", repos = "http://cran.us.r-project.org")
install.packages("tidyverse", repos = "http://cran.us.r-project.org")
install.packages("foreach", repos = "http://cran.us.r-project.org")
install.packages("doParallel", repos = "http://cran.us.r-project.org")
install.packages("parallel", repos = "http://cran.us.r-project.org")


rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)

#macrel.dir="/Users/12705859/bins_clustering_parsing_DFs/macrel_stuff"
#out.dir="/Users/12705859/bins_clustering_parsing_DFs/macrel_stuff/out"
macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"
out.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel/out"

ss1 <- readRDS(file.path(macrel.dir,"ss1.rds"))
ss2 <- readRDS(file.path(macrel.dir,"ss2.rds"))

df1 <- do.call(rbind.data.frame, ss1)
df2 <- do.call(rbind.data.frame, ss2)

df12 <- rbind(df1,df2)

df12$value <- as.numeric(df12$value)
df_sign <- df12 %>% 
  dplyr::filter(value<0.01) %>%
  dplyr::arrange(cluster,speciesB)

head(df_sign)
NROW(df_sign)

#################################################################################
# read correlation data and AMP data 

corr_data <- read.csv(file = file.path(macrel.dir,"corr_data_final"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
corr_data$pig <- as.character(corr_data$pig)

amp_data <- read.csv(file = file.path(macrel.dir,"amp_data_final"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
amp_data$pig <- as.character(amp_data$pig)

melted <- melt(amp_data, id=c("pig","species"), 
     measure.vars=grep("^Clu", colnames(amp_data)))
head(melted)

melted2 <- melted %>%
  dplyr::mutate(cluster = str_replace(variable, "\\."," ")) %>%
  dplyr::select(pig,species,cluster,value)
head(melted2)

#################################################################################

head(corr_data)
head(melted2)


myprint_fun <- function(myDF) {
  
  name <- as.character(myDF$seq[1])
  
  # plot
  pdf(file.path(out.dir,paste0("plot_",name,".pdf")), onefile = TRUE)
  # for row in stats-hits-df, subset and plot : 
  for (row in 1:nrow(myDF)) {
    
    pvalue <- myDF[row,1]
    clu <- myDF[row,2]
    spA <- myDF[row,3]
    spB <- myDF[row,4]
    
    # subset corr
    corr <- subset(corr_data, (speciesA %in% spA))
    corr <- subset(corr, (speciesB %in% spB))
    
    # subset amp
    amp <- subset(melted2, (cluster %in% clu))
    
    # join
    joined <- inner_join(corr,amp, by=c("pig","species")) %>%
      #dplyr::select(speciesA,speciesB,median_corr,pvalue,pig,species,`Cluster 15993`) %>%
      dplyr::mutate_each(funs(replace(., which(is.na(.)), 0))) 
    
    sample_size=NROW(unique(joined$pig))
    
    print(ggplot(joined, aes(x=fct_reorder(pig,as.numeric(joined[,8])),y=median_corr,color=as.factor(joined[,8])))+
            geom_point() +
            ggtitle(label = paste0(joined$speciesA,"\n", joined$cluster), 
                    subtitle = paste0(joined$speciesB,"\n p = ",round(pvalue,8),"\n n = ",sample_size))+
            theme(legend.position="top",
                  legend.title = element_blank(),
                  axis.text.x=element_blank())+
            xlab("pig")+
            ylim(-1,1))
  }
  dev.off()
}


# reduce the number of hits to plots on local machine (or not)
ccc <- df_sign  # df_sign[1:10,] on local

# add sequential numbers
ccc$seq <- seq(1,NROW(ccc))
tail(ccc)

# split into equally sized n groups
multi <- split(ccc, cut(ccc$seq, 60)) # 2 if 2 groups; 60 if 60 groups (these contain 100 each pdf)
NROW(multi)

# apply function (threading): 
cores_to_use <- detectCores()-1
mclapply(multi, myprint_fun, mc.cores = cores_to_use)

