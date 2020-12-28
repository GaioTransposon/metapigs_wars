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



res_out <- read.csv(file = file.path(macrel.dir,"corr_AMPs_stats_out"), sep = ",", 
                    row.names = NULL, header = TRUE, stringsAsFactors = FALSE)

corr_data <- read.csv(file = file.path(macrel.dir,"corr_data_clean"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
corr_data$pig <- as.character(corr_data$pig)

amp_data <- read_csv(file = file.path(macrel.dir,"amp_data.binary_clean"),
                     col_types = cols(pig = col_character()))



# filter out all the hits with pvalues > 0.05 
NROW(res_out)
res_signif <- res_out %>% 
  dplyr::filter(value<0.05)
NROW(res_signif)


# empty df when just one cluster 
res4 <- data.frame(value=double(),
                   cluster=character(),
                   species=character())

res6 <- res_signif   # run res_signif[1:5,] for a quick test
rownames(res6) <- c(seq(1,NROW(res6)))

#pdf(file.path(corr.dir,"amps_to_corrs_sign.pdf"), onefile = TRUE)
for (A in rownames(res6)) {
  
  class(res6)
  A <- as.numeric(A) # dataframes rownames must be taken as numeric
  clu <- res6$cluster[A]
  
  amp_data_sub <- amp_data %>%
    dplyr::select(pig,species,eval(clu)) 
  
  # subsetting of original dataframe based on what is statistically significant (rows of df)
  pp <- full_join(corr_data,amp_data_sub, by=c("pig","species")) %>%
    dplyr::select(speciesA,speciesB,median_corr,pvalue,pig,species,eval(clu)) %>%
    dplyr::mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
    dplyr::filter(., speciesB == as.character(res6$species[A]))
  
  # save some parameters to report on plot
  spB <- unique(pp$speciesB)
  
  zz <- pp %>% 
    dplyr::filter(pp[,7]==1) 
  zz <- as.list(sort(unique(zz$speciesA)))
  
  sp_to_remove <- as.list(unique(pp$speciesB))
  
  check <- subset(pp, (speciesA %in% zz))
  check <- subset(check, (!speciesA %in% sp_to_remove))
  
  head(check)
  
  if ( NROW(check) > 0 & NCOL(check) == 7 & NROW(unique(check[,7]))>1   ) {
    
    
    res0 <- wilcox.test(median_corr ~ check[,7], data=check)
    
    # parse the results:
    value <- res0$p.value
    
    res2 <- data.frame(value,clu,spB)
    res4 <- rbind(res4, res2)
    
  }
  
  else {
    print("no species left to analyze")
  }
}


NROW(res4)

# save the results 
fwrite(x=res4, file = file.path(macrel.dir,"corr_AMPs_stats_out_2nd"))



# pp %>%
#   group_by(`Cluster 15993`) %>%
#   summarise(Unique_speciesA = n_distinct(speciesA),
#             Unique_pigs = n_distinct(pig))



