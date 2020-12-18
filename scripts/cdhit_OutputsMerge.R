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


cdhit.dir = "/shared/homes/12705859/cdhit_work/cd_hit_onAss" # on HPC
#cdhit.dir = "/Users/12705859/Desktop" # on local

clustering_files = list.files(cdhit.dir, pattern="contig_amp_clustering_")

mylist <- lapply(clustering_files, function(x) {
  read.csv(file.path(cdhit.dir,x), sep = ",", 
           row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
  
})

mylist <- lapply(mylist, function(x) { as.data.frame(x) })
mylist <- lapply(mylist, function(x) { x["identity_perc"] <- NULL; x })
mylist <- lapply(mylist, function(x) { x["aa_length"] <- NULL; x })
mylist <- lapply(mylist, function(x) { x["pig"] <- as.character(x$pig); x })

contigs_all_clusters <- Reduce(full_join, mylist)

# move some columns ahead 
contigs_all_clusters <- contigs_all_clusters %>%
  dplyr::select(ORF_number, everything()) %>%
  dplyr::select(contig, everything()) %>%
  dplyr::select(pig, everything()) 


# save all clustering info available for all contigs
fwrite(x = contigs_all_clusters,
       file = file.path(cdhit.dir,"contigs_all_clusters"))


