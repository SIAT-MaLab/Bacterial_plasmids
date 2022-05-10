library(dplyr)
contig_len <- read.csv('contig_len',sep='\t',header=F)
contig_len$breaks <- cut(contig_len$V2,breaks=seq(0,200000,2000))
tmp1 <- contig_len %>% group_by(breaks) %>% sample_n(50,replace=TRUE)