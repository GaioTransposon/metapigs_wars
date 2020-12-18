#upload all libraries
library(base)
library(readr)
library(data.table)
library(dplyr)
library(stringr)
library(utils)
library(reshape)
library(tidyr)
library(reshape2)



corr.dir = "/Users/12705859/Desktop" # on local 
#corr.dir = "/shared/homes/s1/pig_microbiome/fastspar/fastspar_in_per_pig" # on HPC


df <- read.csv(file = file.path(corr.dir,"corr_pvalues_long_concat"), sep = ",", 
               row.names = NULL, header = TRUE, stringsAsFactors = FALSE)

head(df)
# tail(df)
# NROW(df)
# 
# df1 <- df %>% 
#   dplyr::filter(pvalue<0.05)
# 
# NROW(df1)
# 
# df2 <- df1 %>% 
#   dplyr::filter(pvalue<0.01)
# NROW(df2)
# View(df2)
# 
# df2$interaction <- paste0(df2$speciesA,df2$speciesB)
# 
# df3 <- df2 %>%
#   group_by(interaction) %>%
#   dplyr::summarise(mean=mean(median_corr),
#                    sd=sd(median_corr),
#                    n=n()) 





pig <- c(14159, 14159, 14159, 14159, 29951, 29951, 29951, 29951)
contig <- c(265, 327465, 126, 84560, 345867, 8254782, 837, 0239875)
bin <- c('bins.1.fa','bins.10.fa','bins.100.fa','bins.2.fa','bins.1.fa','bins.10.fa','bins.100.fa','bins.2.fa')
species <- c('OTU_A','OTU_B','OTU_C','OTU_D','OTU_A','OTU_B','OTU_C','OTU_D')
cluster100 <- c('clusterX','clusterT','clusterX','clusterW','clusterX','clusterT','clusterY','clusterY')
cluster90 <- c('clusterX','clusterY','clusterZ','clusterW','clusterX','clusterY','clusterX','clusterS')
amp.data <- data.frame(pig, contig, bin, species, cluster100, cluster90)

amp.data$pig_species <- paste0(amp.data$pig,".",amp.data$species)

amp.data_sub <- amp.data %>% 
  dplyr::select(pig_species, cluster100)

si <- unique(amp.data_sub)
amp.data.binary <- dcast(si, formula = pig_species ~ cluster100, fun.aggregate = length)

# split column "pig_species" using _ separator
amp.data.binary <- cSplit(amp.data.binary, "pig_species", ".")

# move last two cols ahead
amp.data.binary <- amp.data.binary %>%
  dplyr::select(pig_species_1, pig_species_2, everything()) 

amp.data.binary <- dplyr::rename(amp.data.binary, pig = pig_species_1)
amp.data.binary <- dplyr::rename(amp.data.binary, species = pig_species_2)

head(amp.data.binary)









pig <- c(14159, 14159, 14159, 14159, 29951, 29951, 29951, 29951,
         14159, 14159, 14159, 14159, 29951, 29951, 29951, 29951,
         14159, 14159, 14159, 14159, 29951, 29951, 29951, 29951,
         14159, 14159, 14159, 14159, 29951, 29951, 29951, 29951)
speciesA <- c('OTU_A','OTU_A','OTU_A','OTU_A','OTU_A','OTU_A','OTU_A','OTU_A',
              'OTU_B','OTU_B','OTU_B','OTU_B','OTU_B','OTU_B','OTU_B','OTU_B',
              'OTU_C','OTU_C','OTU_C','OTU_C','OTU_C','OTU_C','OTU_C','OTU_C',
              'OTU_D','OTU_D','OTU_D','OTU_D','OTU_D','OTU_D','OTU_D','OTU_D')
speciesB <- c('OTU_B','OTU_C','OTU_D','OTU_E','OTU_B','OTU_C','OTU_D','OTU_E',
              'OTU_A','OTU_C','OTU_D','OTU_E','OTU_A','OTU_C','OTU_D','OTU_E',
              'OTU_A','OTU_B','OTU_D','OTU_E','OTU_A','OTU_B','OTU_D','OTU_E',
              'OTU_A','OTU_B','OTU_C','OTU_E','OTU_A','OTU_B','OTU_C','OTU_E')
corr <- c(-0.7, 0.2, 0.3, 0.2, -0.9, 0.1, 0.3, 0.1,
          -0.1, 0.3, 0.3, 0.2, 0.1, 0.1, 0.3, 0.2,
          -0.1, -0.9, 0.3, 0.2, 0.1, 0.9, 0.3, 0.2,
          -0.1, 0.9, 0.3, 0.2, 0.1, 0.9, 0.3, 0.2)
pvalue <- c(0.01, 0.01, 0.01, 0.02, 0.03, 0.02, 0.01, 0.01,
            0.01, 0.01, 0.01, 0.02, 0.03, 0.02, 0.01, 0.01,
            0.03, 0.01, 0.05, 0.02, 0.03, 0.04, 0.01, 0.01,
            0.03, 0.01, 0.05, 0.02, 0.03, 0.04, 0.01, 0.01)

corr.data <- data.frame(pig, speciesA, speciesB, corr, pvalue)

corr.data.edit <- corr.data 
corr.data.edit$species <- corr.data.edit$speciesA
#View(corr.data.edit)




all <- full_join(corr.data.edit,amp.data.binary) %>%
  na.omit(all) # this will not be necessary at a later stage 
head(all)
#View(all)


# split df by speciesB
multiple_DFs <- split( all , f = all$speciesB ,drop = TRUE)
# total number of unique df subsets
NROW(multiple_DFs)

# empty df
res3 <- data.frame(value=double(),
                          cluster=character(),
                          species=character())

# on each subset (based on speciesB) run the following: 
for (single_DF in multiple_DFs) {
  
  #colnames <- colnames(single_DF)
  single <- as.data.frame(single_DF)
  #colnames(single) <- colnames
  
  species <- unique(single$speciesB)
  
  # columns with binary data
  clusters <- single[ , grepl( "cluster" , names( single ) ) ]
  the_rest <- single[ , !grepl( "cluster" , names( single ) ) ]
  
  # remove columns that contain only 0s (otherwise wilcox will give problems cause it won't find two groups to compare)
  clusters <- clusters[colSums(clusters) != 0]
  
  all_clean <- cbind(the_rest,clusters)
  
  # apply independent 2-group Mann-Whitney U Test (designed to work with a continuos vs binary variable):
  res <- lapply(all_clean[ , grepl( "cluster" , names( all_clean ) ) ],
                function(x) wilcox.test(corr ~ x, data=all_clean))
  
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

# save the results 
fwrite(res, x=".....")


res3


# sanity check run test on one and see if it matches with the larger result:
s <- all %>% dplyr::filter(speciesB=="OTU_E")
boxplot(corr ~ clusterX, data = s)
wilcox.test(corr ~ clusterX, data = s)


# what can be done with the output? 

# 1. 
# sort based on pvalue: which clusters look interesting? 
res3 %>% 
  arrange(value)

# 2. 
# filter the original dataframe based on that cluster and species. Visualize. 
all %>% 
  dplyr::select(pig,speciesB,corr,pvalue,clusterX) %>% 
  dplyr::filter(speciesB=="OTU_B") %>% 
  boxplot(corr ~ clusterX, data = .)
