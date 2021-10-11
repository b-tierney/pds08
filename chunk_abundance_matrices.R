#!/usr/bin/Rscript --vanilla

# split pathways and gene family data for parallelization

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

abundance_file = args[[1]]
#abundance_file = 'pathways_baseline.rds'

abfile = readRDS(abundance_file)
abfile = abfile %>% column_to_rownames('Sample_ID') %>% as.data.frame %>% t

chunk = 5000
n = nrow(abfile)
r  = rep(1:ceiling(n/chunk),each=chunk)[1:n]
d = split(as.data.frame(abfile),r)

count = 0
for(file in d){
	count = count + 1
	file = file %>% t %>% as.data.frame %>% rownames_to_column('Sample_ID')
	saveRDS(file,paste(count,'_',abundance_file,sep=''))
}