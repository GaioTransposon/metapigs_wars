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
install.packages("foreach", repos = "http://cran.us.r-project.org")
install.packages("doParallel", repos = "http://cran.us.r-project.org")
install.packages("parallel", repos = "http://cran.us.r-project.org")

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
library(foreach)
library(doParallel)
library(parallel)

#macrel.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/macrel_stuff" # on local 
macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"


corr_data <- read.csv(file = file.path(macrel.dir,"corr_data_clean"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
corr_data$pig <- as.character(corr_data$pig)

amp_data <- read_csv(file = file.path(macrel.dir,"amp_data.binary_clean"),
                     col_types = cols(pig = col_character()))


corr_data2 <- corr_data %>%
  dplyr::filter(speciesA=="Blautia_A sp000285855" | speciesA == "Blautia_A wexlerae" | speciesA == "Faecalibacterium prausnitzii_D"| speciesA == "Faecalibacterium prausnitzii_E") %>%
  dplyr::filter(speciesB=="Agathobacter sp900317585" | speciesB == "Agathobaculum butyriciproducens" | speciesB == "Acutalibacter sp000435395")
NROW(corr_data2)

# split df by both
multiple_DFs <- split(corr_data2, list(corr_data2$speciesB, corr_data2$speciesA), drop = TRUE)
NROW(multiple_DFs)
class(multiple_DFs)



myfun <- function(x) {
  
  colnames <- colnames(x)
  single <- as.data.frame(x)    #single_DF
  colnames(single) <- colnames
  
  speciesA <- as.character(unique(single$speciesA))
  speciesB <- as.character(unique(single$speciesB))
  
  all <- inner_join(single,amp_data, by=c("pig","species"))
  
  all2 <- as.data.frame(all)
  
  clu2 <- all2[,7:ncol(all2)] %>%
    dplyr::mutate_each(funs(replace(., which(is.na(.)), 0)))
  
  all2 <- cbind(all2[,1:6],clu2)
  
  all3 <- all2 %>% 
    dplyr::filter(speciesA!=0) %>% 
    dplyr::filter(speciesB!=0)
  
  evaluated <- all3[,7:ncol(all3)] %>%
    dplyr::mutate(long="long") %>%
    pivot_longer(cols=-long) %>%
    group_by(name) %>%
    dplyr::summarise(sum=sum(value),
                     mean=mean(value)) %>%
    dplyr::filter(sum==0 | mean==1) # if sum is zero it means there s just zeros; if mean is 1 it means there is just 1s
  
  no_pass_list <- as.list(evaluated$name)
  
  all7 <- all3[ , !(names(all3) %in% no_pass_list)]
  
  if ( NROW(all7) > 0 & NCOL(all7) > 7 ) {
    
    res0 <- lapply(all7[ , grepl( "Cluster" , colnames( all7 ) ) ],
                   function(x) wilcox.test(median_corr ~ x, data=all7))
    
    # parse the results:
    res2 <- sapply(res0, function(x) {
      p <- x$p.value
      n <- colnames(p)  
      names(p) <- n
      p
    })
    
    res2 <- melt(res2)
    
    res2$cluster <- rownames(res2)
    rownames(res2) <- NULL
    
    res2$speciesA <- paste0(as.character(speciesA))
    res2$speciesB <- paste0(as.character(speciesB))
    
    return(res2)
    
  }
  
  if (NROW(all7) > 0 & NCOL(all7) == 7) {
    
    res00 <- wilcox.test(median_corr ~ all7[,7], data=all7)
    
    # parse the results:
    value <- res00$p.value
    cluster <- colnames(all7[7])
    speciesA <- paste0(as.character(speciesA))
    speciesB <- paste0(as.character(speciesB))
    
    res22 <- data.frame(value,cluster,speciesA,speciesB)
    
    return(res22)
    
    
  }
  
  if (NROW(all7) > 0 & NCOL(all7) == 6) {
    print(" amp absent in this species ")
  }
  
  else {
    print("something went wrong")
  }
}

cores_to_use <- detectCores()-1
ss <- mclapply(multiple_DFs, myfun, mc.cores = cores_to_use)

saveRDS(ss, file = file.path(macrel.dir,"ss.rds"))

# system.time(
#   ss <- lapply(multiple_DFs[1:2], myfun1)
# ) # HPC; elapsed 108.664 
# 
# system.time(
#   ss <- mclapply(multiple_DFs[1:2], myfun1, mc.cores = 2)
# ) # HPC; elapsed 61.297 
# 
# system.time(
#   ss <- mclapply(multiple_DFs[1:5], myfun1, mc.cores = 10)
# ) # HPC; elapsed 61.650 
# 
# system.time(
#   ss <- mclapply(multiple_DFs[1:10], myfun1, mc.cores = 10)
# ) # HPC; elapsed 67.989 


#rrr <- readRDS("/Users/12705859/Desktop/ss.rds")
rrr <- readRDS("/Users/12705859/Desktop/ss.rds")
library(dplyr)
library(ggplot2)
library(tidyverse) 
df <- ldply(rrr, data.frame)
df <- df %>% 
  dplyr::filter(value<0.01)
df

full_join(corr_data,amp_data.binary3, by=c("pig","species")) %>%
  dplyr::select(speciesA,speciesB,median_corr,pvalue,pig,species,`Cluster 17298`) %>%
  dplyr::mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
  dplyr::filter(., speciesB == "Agathobaculum butyriciproducens") %>%
  dplyr::filter(., speciesA == "Blautia_A wexlerae") %>%
  ggplot(., aes(x=fct_reorder(pig,as.numeric(.[,7])),y=median_corr,color=.[,7]))+
  geom_point() +
  facet_grid(~speciesA, scales = "free")+
  theme(legend.position="top",
        axis.text.x=element_blank())+
  xlab("pig")

full_join(corr_data,amp_data.binary3, by=c("pig","species")) %>%
  dplyr::select(speciesA,speciesB,median_corr,pvalue,pig,species,`Cluster 15993`) %>%
  dplyr::mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
  dplyr::filter(., speciesB == "Agathobacter sp900317585") %>%
  dplyr::filter(., speciesA == "Faecalibacterium prausnitzii_D") %>%
  ggplot(., aes(x=fct_reorder(pig,as.numeric(.[,7])),y=median_corr,color=.[,7]))+
  geom_point() +
  facet_grid(~speciesA, scales = "free")+
  theme(legend.position="top",
        axis.text.x=element_blank())+
  xlab("pig")
