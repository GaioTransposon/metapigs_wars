
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



# subset both based on species and clusters present in my significant-hits output
ss1 <- readRDS("ss1.rds")
ss2 <- readRDS("ss2.rds")

df1 <- do.call(rbind.data.frame, ss1)
df2 <- do.call(rbind.data.frame, ss2)


df1$value <- as.numeric(df1$value)
View(df1_sign)
df1_sign <- df1 %>% 
  dplyr::filter(value<0.01) %>%
  dplyr::arrange(cluster,speciesA)


df2$value <- as.numeric(df2$value)
df2_sign <- df2 %>% 
  dplyr::filter(value<0.01) %>%
  dplyr::arrange(value)


df <- rbind(df1_sign,df2_sign)

head(df)
tail(df)
NROW(df)

NROW(unique(df$speciesA))
NROW(unique(df$speciesB))
NROW(unique(df$cluster))

spA <- as.list(unique(df$speciesA))
spB <- as.list(unique(df$speciesB))
clu <- as.list(unique(df$cluster))

NROW(corr_data3)
corr_data_final <- subset(corr_data3, (speciesA %in% spA))
corr_data_final <- subset(corr_data_final, (speciesB %in% spB))
NROW(corr_data_final)


NCOL(amp_data)
amp_data_final <- amp_data[names(amp_data) %in% clu]
amp_data_final <- cbind(amp_data[,1:2],amp_data_final)
NCOL(amp_data_final)


# save the results 
fwrite(x=corr_data_final, file = file.path(macrel.dir,"corr_data_final"))
# save the results 
fwrite(x=amp_data_final, file = file.path(macrel.dir,"amp_data_final"))
