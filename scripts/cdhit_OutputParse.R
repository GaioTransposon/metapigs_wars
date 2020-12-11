# runs from the HPC 
# language: R 

#This script requires the following packages:
install.packages("base", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org", dependencies = TRUE)
install.packages("dplyr", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("utils", repos = "http://cran.us.r-project.org")
install.packages("splitstackshape", repos = "http://cran.us.r-project.org")

#upload all libraries
library(base)
library(data.table)
library(dplyr)
library(stringr)
library(utils)
library(splitstackshape)


#cdhit.dir = "/shared/homes/12705859/cdhit_work/cd_hit_onAss" # on HPC
cdhit.dir = "/Users/12705859/Desktop" # on local


raw_clustering_files = list.files(cdhit.dir, pattern=".clstr")
for (clstr in raw_clustering_files) {
  
  clstr1 <- read.csv(file.path(cdhit.dir,clstr), sep = "\t", 
                     row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
  
  
  clstr2 <- clstr1
  n = nrow(clstr1)
  x = 0
  numbers_only <- function(x) !grepl("\\D", x)
  for (row in c(1:n)) {
    if (numbers_only(clstr2[row,1]) == TRUE) {
      clstr2[row,1] <- x}
    else {NULL}
    x <- clstr2[row,1]
  }
  
  clstr.sums <- data.frame(dplyr::count(clstr2,V1))
  
  switch <- clstr.sums[1,2]
  clstr3 <- cbind(clstr2[1], clstr1)
  
  clstr3[c((switch-5):(switch+5)),]
  
  clstr4 <- clstr2[-which(clstr2$V2 == ""), ]
  clstr4[c(1:5,(switch-5):(switch+5)),]
  
  clstr5 <- clstr4
  clstr5[] <- lapply(clstr5, gsub, pattern='>', replacement='')
  clstr5.2 <- data.frame(str_split_fixed(clstr5$V2, "aa, ", 2))
  clstr5.3 <- data.frame(str_split_fixed(clstr5.2$X2, "... ", 2))
  clstr6 <- cbind(clstr5[1],clstr5.2[1],clstr5.3[1:2])
  colnames(clstr6) <- c("cluster","aa","gene","stat")
  
  
  # split column "gene" 
  clstr7 <- cSplit(clstr6, "gene", "_")
  
  # rename cols
  colnames(clstr7) <- c("cluster","aa_length","identity_perc","pig","bin","contig")
  
  clstr7$identity_perc <- gsub( "at", "", clstr7$identity_perc)
  clstr7$identity_perc <- gsub( " ", "", clstr7$identity_perc)
  clstr7$identity_perc <- gsub( "%", "", clstr7$identity_perc)
  
  # get minimum percentage of identity across df (to save as affix to output filename)
  sub <- clstr7 %>% 
    dplyr::filter(!identity_perc=="*")
  sub$identity_perc <- as.numeric(sub$identity_perc)
  percentage <- min(sub$identity_perc)
  
  perc_round <- round(min(sub$identity_perc))
  colnames(clstr7)[colnames(clstr7)=="cluster"] <- paste0("cluster_",perc_round)
  
  # save data 
  fwrite(x = clstr7,
         file = file.path(cdhit.dir,paste0("contig_amp_clustering_",percentage,"min_percent")))

}




