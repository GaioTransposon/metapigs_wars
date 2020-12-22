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
tail(df)
NROW(df)

# ## produce some stats about the correlations:
# df1 <- df %>%
#   dplyr::filter(pvalue<0.05)
# NROW(df1)
# 
# df2 <- df1 %>%
#   dplyr::filter(pvalue<=0.01)
# NROW(df2)
# 
# df1$interaction <- paste0(df1$speciesA,df1$speciesB)
# df3 <- df1 %>%
#   group_by(interaction) %>%
#   dplyr::summarise(mean=mean(median_corr),
#                    sd=sd(median_corr),
#                    n=n()) %>%
#   dplyr::arrange(desc(n))
# head(df3)
# NROW(df3)

#######

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
class(amp.data_sub)
si <- unique(amp.data_sub)
amp.data.binary <- dcast(si, formula = pig_species ~ cluster100, fun.aggregate = length)


# remove columns which colSums==1
clusters1<- amp.data.binary[,2:ncol(amp.data.binary)]
NCOL(clusters1)
clusters1 <- clusters1[,-(which(colSums(clusters1)<2))]
NCOL(clusters1)
colSums(clusters1)

amp.data.binary2 <- cbind(amp.data.binary$pig_species,clusters1)
names(amp.data.binary2)[names(amp.data.binary2) == 'amp.data.binary$pig_species'] <- 'pig_species'

# split column 
amp.data.binary2 <- cSplit(amp.data.binary2, "pig_species", ".")

# move last two cols ahead
amp.data.binary3 <- amp.data.binary2 %>%
  dplyr::select(pig_species_1, pig_species_2, everything()) 
names(amp.data.binary3)[names(amp.data.binary3) == 'pig_species_1'] <- 'pig'
names(amp.data.binary3)[names(amp.data.binary3) == 'pig_species_2'] <- 'species'


head(amp.data.binary3)









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
  clusterss <- single[ , grepl( "cluster" , names( single ) ) ]
  the_rest <- single[ , !grepl( "cluster" , names( single ) ) ]
  
  # remove columns that contain only 0s (otherwise wilcox will give problems cause it won't find two groups to compare)
  clusterss <- clusterss[colSums(clusterss) != 0]
  
  all_clean <- cbind(the_rest,clusterss)
  
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
