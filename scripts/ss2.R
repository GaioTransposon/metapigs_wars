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

# remove rows where speciesA and speciesB are identical 
NROW(corr_data)
i<-which(!(corr_data$speciesA==corr_data$speciesB))
corr_data<-corr_data[i, ]
NROW(corr_data)

# corr_data2 <- corr_data %>%
#   dplyr::filter(speciesA=="Blautia_A sp000285855" | speciesA == "Blautia_A wexlerae" | speciesA == "Faecalibacterium prausnitzii_D"| speciesA == "Faecalibacterium prausnitzii_E") %>%
#   dplyr::filter(speciesB=="Agathobacter sp900317585" | speciesB == "Agathobaculum butyriciproducens" | speciesB == "Acutalibacter sp000435395")
# NROW(corr_data2)

corr_data2 <- corr_data

corr_data2$speciesAspeciesB <- paste0(corr_data2$speciesA,".",corr_data2$speciesB)

keep <- corr_data2 %>% 
  group_by(speciesAspeciesB) %>%
  dplyr::summarise(sum=n()) %>% 
  dplyr::filter(sum>4)

NROW(corr_data2)
corr_data3 <- subset(corr_data2, speciesAspeciesB %in% keep$speciesAspeciesB )
NROW(corr_data3)

corr_data3$speciesAspeciesB <- NULL

# split df by both
multiple_DFs <- split(corr_data3, list(corr_data3$speciesB, corr_data3$speciesA), drop = TRUE)
NROW(multiple_DFs)

#multiple_DFs1 <- multiple_DFs[1:129102]
multiple_DFs2 <- multiple_DFs[129103:NROW(multiple_DFs)]


myfun <- function(x) {
  
  colnames <- colnames(x)
  single <- as.data.frame(x)    #single_DF
  colnames(single) <- colnames
  
  speciesA <- as.character(unique(single$speciesA))
  speciesB <- as.character(unique(single$speciesB))
  
  # full_join produces output on HPC, testing 10 items of the list
  all <- full_join(single,amp_data, by=c("pig","species"))
  
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
ss2 <- mclapply(multiple_DFs2, myfun, mc.cores = cores_to_use)

saveRDS(ss2, file = file.path(macrel.dir,"ss2.rds"))
