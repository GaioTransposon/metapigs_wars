
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

#library(foreach)
#library(doParallel)
#library(purrr)
library(parallel)


#macrel.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/macrel_stuff" # on local 
macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"


corr_data <- read.csv(file = file.path(macrel.dir,"corr_data_clean"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
corr_data$pig <- as.character(corr_data$pig)

amp_data <- read_csv(file = file.path(macrel.dir,"amp_data.binary_clean"),
                     col_types = cols(pig = col_character()))
amp_data = data.table(amp_data, key=c("pig","species"))
NROW(amp_data)
NCOL(amp_data)

# remove rows where speciesA == speciesB
NROW(corr_data)
NROW(unique(corr_data$speciesB))
corr_data2 <- subset(corr_data, speciesA != speciesB)
NROW(corr_data2)
# keep species present in the amp dataset
corr_data3 <- corr_data2[corr_data2$speciesA %in% amp_data$species,]
NROW(corr_data3)


# split df by speciesA and speciesB
multiple_DFs_double <- split(corr_data3,list(corr_data3$speciesB,corr_data3$speciesA), drop = TRUE)

# total number of unique df subsets
NROW(multiple_DFs_double)
test <- multiple_DFs_double    # multiple_DFs_double[1:50]
NROW(test)

my_join_fun <- function(single_DF) {
  
  colnames <- colnames(single_DF)
  single <- as.data.frame(single_DF)
  colnames(single) <- colnames
  
  speciesA <- as.character(unique(single$speciesA))
  speciesB <- as.character(unique(single$speciesB))
  
  single$species <- single$speciesA
  
  all0 <- inner_join(single, amp_data, by =c("pig","species") )
  
  all1 <- list(all0)
  return(all1)

}

NROW(test)
df_list <- NULL
#df_list <- sapply(test, my_join_fun)
#df_list <- mclapply(test, my_join_fun, mc.cores=1)
df_list <- mclapply(test, my_join_fun, mc.cores=40)
NROW(df_list)


to_save <- do.call("rbind", lapply(df_list, as.data.frame)) 


# 
# 
# 
# df_list2 <- df_list[sapply(df_list, function(x) dim(x)[1]) > 1]
# NROW(df_list2)
# 
# 
# 
# my_NAtoZERO_fun <- function(df) {
#   
#   #df <- df_list2$`Solobacterium timonensis.Butyricicoccus_A sp002395695`
#   
#   df <- df %>%
#     dplyr::mutate_each(funs(replace(., which(is.na(.)), 0)))
#   
#   clusters_to_remove <- melt(df[,7:ncol(df)]) %>% 
#     group_by(variable) %>%
#     dplyr::summarise(Unique_Elements = n_distinct(value)) %>% 
#     dplyr::filter(Unique_Elements==1) 
# 
#   all7 <- df[ , !(names(df) %in% clusters_to_remove$variable)]
#   
#   return(all7)
#   
# }
# 
# NROW(df_list2)
# df_list3 <- NULL
# df_list3 <- sapply(df_list2, my_NAtoZERO_fun)
# NROW(df_list3)
# 
# df_list4 <- Filter(function(x) nrow(x) > 0, df_list3)
# df_list5 <- Filter(function(x) ncol(x) > 6, df_list4)
# NROW(df_list5)
# 
# # once we get a df with some data (after the above filtering we set up wilcox)
# 
# ######
# lapply(df_list6[ , grepl( "Cluster" , colnames( df_list6 ) ) ],
#        function(x) wilcox.test(median_corr ~ x, data=df_list6))
# 
# 
# 
# 
# 
# 
# 
#   #print(colnames(all))
#   print(NROW(all))
#   
#   if (NROW(all) > 4) {
#     
#     clu <- all[,7:ncol(all)] %>%
#       dplyr::mutate_each(funs(replace(., which(is.na(.)), 0)))
#     
#     all3 <- cbind(all[,1:6],clu)
#     
#     to_evaluate <- all3[,7:ncol(all3)] %>%
#       dplyr::summarise_all(funs(sum))
#     to_evaluate$long <- paste0("long")
#     remove_if_tot_sum_less_than_five <- to_evaluate %>% 
#       pivot_longer(cols=-long) %>%
#       dplyr::filter(value<5)
#     no_pass_list <- as.list(remove_if_tot_sum_less_than_five$name)
#     
#     all7 <- all3[ , !(names(all3) %in% no_pass_list)]
#     
#     if ( NROW(all7) > 0 & NCOL(all7) > 7 ) {
#       
#       
#       res0 <- lapply(all7[ , grepl( "Cluster" , colnames( all7 ) ) ],
#                      function(x) wilcox.test(median_corr ~ x, data=all7))
#       
#       # parse the results:
#       res2 <- sapply(res0, function(x) {
#         p <- x$p.value
#         n <- colnames(p)  
#         names(p) <- n
#         p
#       })
#       
#       res2 <- melt(res2)
#       
#       res2$cluster <- rownames(res2)
#       rownames(res2) <- NULL
#       
#       res2$speciesA <- paste0(as.character(speciesA))
#       res2$speciesB <- paste0(as.character(speciesB))
#       
#       # rbind the results
#       res3 <- rbind(res3,res2)
#       
#     }
#     
#     if (NROW(all7) > 0 & NCOL(all7) == 7) {
#       
#       res00 <- wilcox.test(median_corr ~ all7[,7], data=all7)
#       
#       # parse the results:
#       value <- res00$p.value
#       cluster <- colnames(all7[7])
#       speciesA <- paste0(as.character(speciesA))
#       speciesb <- paste0(as.character(speciesB))
#       
#       res22 <- data.frame(value,cluster,speciesA, speciesB)
#       
#       res4 <- rbind(res4, res22)
#       
#     }
#     
#     else {
#       print("contig not associating to this species")
#     }
#   }
#   
#   else {
#     print("no amp data for this species")
#   }
#   
# }
# 
# 
#   
#   
# 
#   
#   
# 
# res5 <- rbind(res3,res4) %>%
#   dplyr::arrange(value)
# head(res5)
# 
# NROW(res5)
# 
# 
# NROW(unique(mtcars$cyl))
# 
# 
# z <- split(mtcars,list(mtcars$cyl,mtcars$gear))
# NROW(z)
