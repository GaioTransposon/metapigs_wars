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


#corr.dir = "/Users/12705859/Desktop" # on local 
corr.dir = "/shared/homes/s1/pig_microbiome/fastspar/fastspar_in_per_pig" # on HPC

#macrel.dir = "/Users/12705859/Desktop" # on local 
macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"


corr_data <- read.csv(file = file.path(corr.dir,"corr_pvalues_long_concat"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)

corr_data.edit <- corr_data 

sort(unique(corr_data.edit$speciesB))
# remove NA species
NROW(corr_data.edit)
corr_data.edit2 <- corr_data.edit[!grepl("NA__",corr_data.edit$speciesA),]
corr_data.edit2 <- corr_data.edit2[!grepl("NA__",corr_data.edit2$speciesB),]
NROW(corr_data.edit2)
corr_data.edit2$pig <- as.character(corr_data.edit2$pig)

corr_data.edit2$speciesA <- substr(corr_data.edit2$speciesA, 1, regexpr("\\__", corr_data.edit2$speciesA)-1)
corr_data.edit2$speciesB <- substr(corr_data.edit2$speciesB, 1, regexpr("\\__", corr_data.edit2$speciesB)-1)
corr_data.edit2$speciesA <- sub("^X", "", corr_data.edit2$speciesA)
corr_data.edit2$speciesB <- sub("^X", "", corr_data.edit2$speciesB)
corr_data.edit2$speciesA <- gsub("\\-", " ", corr_data.edit2$speciesA)
corr_data.edit2$speciesB <- gsub("\\.", " ", corr_data.edit2$speciesB)
head(corr_data.edit2)
tail(corr_data.edit2)

corr_data.edit2$species <- corr_data.edit2$speciesA

a <- sort(unique(corr_data.edit2$speciesA))
b <- sort(unique(corr_data.edit2$speciesB))
a <- as.data.frame(a)
b <- as.data.frame(b)
# if TRUE, same now match
NROW(a)==NROW(which(a==b))
###



###
merge_pig_contig_bin_species_clusters_out <- read_csv(file = file.path(macrel.dir,"merge_pig_contig_bin_species_clusters_out"),
                                                      col_types = cols(pig = col_character()))

amp_data <- as.data.frame(merge_pig_contig_bin_species_clusters_out)


### use the amp data to find which species are not present in the corr data, and remove them 
d <- sort(unique(amp_data$species))
d <- as.data.frame(d)
d$d <- gsub("\\-", " ", d$d)
head(d)
head(a)
colnames(d) <- "a"
# all the species from the corr dataframe that are not present in the AMP table: 
not_present <- as.list(anti_join(a, d, by = c("a")))

# remove them from the corr dataframe: 
corr_part1 <- subset(corr_data.edit2, !(speciesB %in% not_present$a))
corr_part2 <- subset(corr_part1, !(speciesA %in% not_present$a))
head(corr_data.edit2)
NROW(corr_data.edit2)
NROW(corr_part1)
NROW(corr_part2)
###

# continue with the amp data 
amp_data$species <- gsub("\\-", " ", amp_data$species)



test <- amp_data %>% 
  dplyr::select(pig, species, cluster_98) %>%
  group_by(pig,cluster_98) %>% 
  tally() %>%
  dplyr::arrange(desc(n))

test2 <- amp_data %>% 
  dplyr::select(pig, species, cluster_98) %>%
  group_by(cluster_98) %>% 
  tally() %>%
  dplyr::arrange(desc(n))

#############################################
#############################################
## IMPORTANT TO RUN ON LOCAL
#############################################
#############################################
# amp_data <- amp_data %>% 
#   dplyr::filter(cluster_98=="Cluster 22071"|
#                   cluster_98=="Cluster 15993"|
#                   cluster_98=="Cluster 11017"|
#                   cluster_98=="Cluster 25291"|
#                   cluster_98=="Cluster 19756")
#############################################
#############################################

amp_data$pig_species <- paste0(amp_data$pig,".",amp_data$species)
amp_data_sub <- amp_data %>% 
  dplyr::select(pig_species, cluster_98) 


si <- unique(amp_data_sub)
NROW(si)
# remove species that contain NA (these contigs don't associate to any species so it's worthless taking them on to look at species corr)
si <- si[!grepl(".NA",si$pig_species),]
NROW(si)


amp_data.binary <- reshape2::dcast(si, formula = pig_species ~ cluster_98, fun.aggregate = length)

# split column 
amp_data.binary2 <- cSplit(amp_data.binary, "pig_species", ".")
# move last two cols ahead
amp_data.binary3 <- amp_data.binary2 %>%
  dplyr::select(pig_species_1, pig_species_2, everything()) 
names(amp_data.binary3)[names(amp_data.binary3) == 'pig_species_1'] <- 'pig'
names(amp_data.binary3)[names(amp_data.binary3) == 'pig_species_2'] <- 'species'
amp_data.binary3$pig <- as.character(amp_data.binary3$pig)
##


dim(corr_part2)
dim(amp_data.binary3)
head(corr_part2)
head(amp_data.binary3)

# apply stats on subset:
# split df by speciesB
multiple_DFs <- split( corr_part2 , f = corr_part2$speciesB, drop = TRUE)
# total number of unique df subsets
NROW(multiple_DFs)
names(multiple_DFs)

# empty df
res3 <- data.frame(value=double(),
                   cluster=character(),
                   species=character())



# on each subset (based on speciesB) run the following: 
for (single_DF in multiple_DFs) {
  
  colnames <- colnames(single_DF)
  single <- as.data.frame(single_DF)    #single_DF
  colnames(single) <- colnames
  
  species <- as.character(unique(single$speciesB))
  
  
  all <- full_join(single,amp_data.binary3, by=c("pig","species"))
  
  
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
  
  
  ## this below was not working on HPC 
  #all5 <- all4[c(rep(TRUE)), colSums(all4[7L:ncol(all4)]) > 0L]
  
  
  if ( NROW(all7) > 0 & NCOL(all7) > 6 ) {
    
    res <- lapply(all7[ , grepl( "Cluster" , colnames( all7 ) ) ],
                  function(x) wilcox.test(median_corr ~ x, data=all7))
    
    # this below it's another equally valid stat test. The output needs to be parsed differently. 
    # res <- lapply(all_clean[ , grepl( "cluster" , names( all_clean ) ) ],
    #               function(x) glm(x ~ corr, data=all_clean))
    
    # parse the results:
    res2 <- sapply(res, function(x) {
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
  
  else {
    print(" contig not associating to this species ")
  }
  
}

res4 <- res3 %>%
  dplyr::arrange(value)
head(res4)

# save the results 
fwrite(x=res4, file = file.path(macrel.dir,"corr_AMPs_stats_out"))



## vis stuff
# #library(ggplot2)
# df1 <- full_join(corr_part2,amp_data.binary3, by=c("pig","species")) %>% 
#   dplyr::mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
#   dplyr::filter(speciesB=="Acetatifactor sp900066365") %>%
#   dplyr::select(speciesA,speciesB,median_corr,pvalue,pig,species,`Cluster 15993`)
# df1 %>%
#   ggplot(., aes(x=as.character(`Cluster 15993`), y=median_corr))+
#   geom_boxplot()
# hist(df1$median_corr)
# 
# df1 <- full_join(corr_part2,amp_data.binary3, by=c("pig","species")) %>% 
#   dplyr::mutate_each(funs(replace(., which(is.na(.)), 0))) %>%
#   dplyr::filter(speciesB=="14 2 sp001940225") %>%
#   dplyr::select(speciesA,speciesB,median_corr,pvalue,pig,species,`Cluster 11017`)
# df1 %>%
#   ggplot(., aes(x=as.character(`Cluster 11017`), y=median_corr))+
#   geom_boxplot()
# hist(df1$median_corr)




## evt to run on all without subsetting (to work on)
# 
# single <- as.data.frame(corr_data.edit)    #single_DF
# 
# all <- full_join(single,c) 
# 
# # ...
# all2 <- all %>% 
#   dplyr::mutate_each(funs(replace(., which(is.na(.)), 0)))
# 
# # remove columns that contain only 0s (otherwise wilcox will give problems cause it won't find two groups to compare)
# clusters <- all2[,colSums(all2[,7:ncol(all2)]) > 2]
# 
# cc <- cbind(all2[,1:6],clusters)
# 
# res <- lapply(cc[ , grepl( "Cluster" , names( cc ) ) ],
#               function(x) wilcox.test(median_corr ~ x, data=cc))
# 
# # parse the results:
# res2 <- sapply(res, function(x) {
#   p <- x$p.value
#   n <- colnames(p)  
#   names(p) <- n
#   p
# })
# 
# res2 <- melt(res2)
# res2$cluster <- rownames(res2)
# rownames(res2) <- NULL
# res2$species <- paste0(as.character(species))
# 
# res3 <- res2 %>%
#   dplyr::arrange(value)
# 
# # save the results 
# fwrite(x=res3, file = file.path(macrel.dir,"corr_AMPs_stats_out98_clu"))


