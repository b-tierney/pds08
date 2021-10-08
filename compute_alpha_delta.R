#!/usr/bin/Rscript --vanilla

# compute alpha delta

#setwd('/Volumes/GoogleDrive-109433674545306273960/My\ Drive/pds08')
setwd('~/Desktop/gdrive/pds08')

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

base = args[[1]]
end = args[[2]]

basename = base %>% strsplit('_',base) %>% map_chr(1)

basedata = readRDS(base) 
endata = readRDS(end) 

basedata = basedata %>% filter(Sample_ID %in% endata$Sample_ID) %>% arrange(Sample_ID) %>% column_to_rownames('Sample_ID')
endata = endata %>% filter(Sample_ID %in% rownames(basedata)) %>% arrange(Sample_ID) %>% column_to_rownames('Sample_ID')

delta = basedata - endata 
delta = delta %>% rownames_to_column('Sample_ID')

saveRDS(delta,paste(basename,'_delta_diversity.rds',sep=''))