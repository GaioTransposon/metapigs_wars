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
library(ggplot2)
library(EnvStats)
library(ggpubr)
library(patchwork)
library(tidyverse) # important for fct function in last plot (bottom_plot)


macrel.dir = "/Users/12705859/Desktop/bins_clustering_parsing_DFs/macrel_stuff" # on local 
#macrel.dir="/shared/homes/12705859/cdhit_work/cd_hit_onMacrel"


corr_data <- read.csv(file = file.path(macrel.dir,"corr_data_clean"), sep = ",", 
                      row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
corr_data$pig <- as.character(corr_data$pig)

amp_data <- read_csv(file = file.path(macrel.dir,"amp_data.binary_clean"),
                     col_types = cols(pig = col_character()))

res_out <- read.csv(file = file.path(macrel.dir,"corr_AMPs_stats_out_2nd"), sep = ",", 
                    row.names = NULL, header = TRUE, stringsAsFactors = FALSE)

names(res_out)[names(res_out) == 'spB'] <- 'species'

# filter out all the hits with pvalues > 0.05 
NROW(res_out)
res_signif <- res_out %>% 
  dplyr::filter(value<0.05)
NROW(res_signif)

# pp %>%
#   group_by(`Cluster 15993`) %>%
#   summarise(Unique_speciesA = n_distinct(speciesA),
#             Unique_pigs = n_distinct(pig))


# empty df when just one cluster 
res4 <- data.frame(value=double(),
                   cluster=character(),
                   species=character())

res6 <- res_signif   # run res_signif[1:5,] for a quick test
rownames(res6) <- c(seq(1,NROW(res6)))

pdf(file.path(macrel.dir,"amps_to_corrs_sign.pdf"), onefile = TRUE)
for (A in rownames(res6)) {
  
  class(res6)
  A <- as.numeric(A) # dataframes rownames must be taken as numeric
  clu <- res6$clu[A]
  
  # subsetting of original dataframe based on what is statistically significant (rows of df)
  pp <- full_join(corr_data,amp_data, by=c("pig","species")) %>%
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
  
  # # build plot
  # plot <- ggplot(pp, aes(x=as.character(pp[,7]), y=median_corr))+
  #   geom_boxplot()+
  #   #geom_point(alpha=0.3) +
  #   stat_n_text() +
  #   theme(axis.title.x = element_text(size=15),
  #         axis.title.y = element_text(size=15))+
  #   ggtitle(spB,
  #           subtitle= as.character(eval(clu))) +
  #   xlab(NULL)
  # 
  # dens <- ggplot(pp, aes(x = median_corr, fill = as.factor(pp[,7]))) +
  #   geom_density(alpha = 0.4) +
  #   xlab("")+
  #   theme(legend.position = "right",
  #         legend.title = element_blank())+
  #   coord_flip()

  last <- ggplot(check, aes(x=as.character(check[,7]), y=median_corr))+
    geom_boxplot()+
    stat_n_text() +
    theme(axis.title.y = element_text(size=10))+
    xlab(NULL) +
    ggtitle(spB,
            subtitle= as.character(eval(clu)))
  
  check[,7] <- as.factor(check[,7])

  NROW(unique(check$speciesA))
  
  bottom_plot <- check %>%
    ggplot(., aes(x=fct_reorder(pig,as.numeric(check[,7])),y=median_corr,color=check[,7]))+
    geom_point() +
    facet_grid(~speciesA, scales = "free")+
    theme(legend.position="none",
          axis.text.x=element_blank())+
    xlab("pig")
  
  final <- last + bottom_plot +
    plot_layout(ncol = 2,
                nrow = 1,
                widths = c(3, 4))
  
  print(final)
}
dev.off()





