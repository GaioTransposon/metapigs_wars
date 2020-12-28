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


macrel.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/macrel_stuff" # on local 
#macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"


corr_data <- read.csv(file = file.path(macrel.dir,"corr_data_clean"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
corr_data$pig <- as.character(corr_data$pig)

amp_data <- read_csv(file = file.path(macrel.dir,"amp_data.binary_clean"),
                                                      col_types = cols(pig = col_character()))

# apply stats on subset:
# split df by speciesB
multiple_DFs <- split( corr_data , f = corr_data$speciesB, drop = TRUE)

# total number of unique df subsets
NROW(multiple_DFs)

# empty df when more than 1 cluster
res3 <- data.frame(value=double(),
                   cluster=character(),
                   species=character())

# empty df when just one cluster 
res4 <- data.frame(value=double(),
                   cluster=character(),
                   species=character())

# on each subset (based on speciesB) run the following: 
for (single_DF in multiple_DFs[1:300]) {
  
  colnames <- colnames(single_DF)
  single <- as.data.frame(single_DF)    #single_DF
  colnames(single) <- colnames
  
  species <- as.character(unique(single$speciesB))
  
  print(species)
  
  all <- full_join(single,amp_data, by=c("pig","species"))
  
  
  clu <- all[,7:ncol(all)] %>%
    dplyr::mutate_each(funs(replace(., which(is.na(.)), 0)))
  
  all2 <- cbind(all[,1:6],clu)
  
  all3 <- all2 %>% 
    dplyr::filter(speciesB!=0)
  
  to_evaluate <- all3[,7:ncol(all3)] %>%
    dplyr::summarise_all(funs(sum))
  to_evaluate$long <- paste0("long")
  
  evaluated <- to_evaluate %>% 
    pivot_longer(cols=-long) %>%
    dplyr::filter(value<5)
  
  no_pass_list <- as.list(evaluated$name)
  
  all7 <- all3[ , !(names(all3) %in% no_pass_list)]
  
  print(names(all3)[1:10])
  print(names(all7)[1:10])
  ## this below was not working on HPC 
  #all5 <- all4[c(rep(TRUE)), colSums(all4[7L:ncol(all4)]) > 0L]
  
  
  if ( NROW(all7) > 0 & NCOL(all7) > 7 ) {
    
    
    res0 <- lapply(all7[ , grepl( "Cluster" , colnames( all7 ) ) ],
                  function(x) wilcox.test(median_corr ~ x, data=all7))
    
    # # this below it's another equally valid stat test. The output needs to be parsed differently. 
    # res <- lapply(all7[ , grepl( "cluster" , names( all7 ) ) ],
    #               function(x) glm(median_corr ~ x, data=all7))
    # 
    # glm(median_corr ~ all7$`Cluster 11017`, data=all7)
    
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
    
    res2$species <- paste0(as.character(species))
    
    # rbind the results
    res3 <- rbind(res3,res2)
    
  }
  
  if (NROW(all7) > 0 & NCOL(all7) == 7) {
    
    res00 <- wilcox.test(median_corr ~ all7[,7], data=all7)
    
    # parse the results:
    value <- res00$p.value
    cluster <- colnames(all7[7])
    species <- paste0(as.character(species))
    
    res22 <- data.frame(value,cluster,species)
    
    res4 <- rbind(res4, res22)
    
  }

  else {
    print(" contig not associating to this species ")
  }
  
}

res5 <- rbind(res3,res4) %>%
  dplyr::arrange(value)
head(res5)
NROW(res5)

# save the results 
fwrite(x=res5, file = file.path(macrel.dir,"corr_AMPs_stats_out"))

