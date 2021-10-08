#!/usr/bin/Rscript --vanilla

# compute alpha and beta diversities

library(tidyverse)
library(vegan)

setwd('/Volumes/GoogleDrive-109433674545306273960/My\ Drive/pds08')

args <- commandArgs(trailingOnly = TRUE)

dataloc = args[[1]]
basename = gsub('.rds','',dataloc)
data = readRDS(dataloc)

data_t = data %>% column_to_rownames('Sample_ID') %>% t
shannon = map(colnames(data_t), function(x) diversity(data_t[,x],index='shannon')) %>% unlist %>% unname
simpson = map(colnames(data_t), function(x) diversity(data_t[,x],index='simpson')) %>% unlist %>% unname
data_t_species = data_t[grepl('s__',rownames(data_t)),]
richness = map(colnames(data_t_species), function(x) length(which(data_t_species[,x]!=0))) %>% unlist %>% unname

divs = data.frame(shannon,simpson,richness)
rownames(divs) = colnames(data_t)
divs = divs %>% rownames_to_column('Sample_ID')

saveRDS(divs,paste(basename,'_diversity.rds',sep=''))