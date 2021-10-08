#!/usr/bin/Rscript --vanilla

# compute associations with phenotypes

library(tidyverse)
library(broom)

#setwd('/Volumes/GoogleDrive-109433674545306273960/My\ Drive/pds08')
setwd('~/Desktop/gdrive/pds08')

args <- commandArgs(trailingOnly = TRUE)

metadata_file = args[[1]]
abundance_data_file = args[[2]]

#metadata_file='pds08_metadata.rds'
#abundance_data_file = 'metaphlan_endpoint_diversity.rds'

metadata = readRDS(metadata_file) 
abundance_data = readRDS(abundance_data_file)

#abundance_data = abundance_data %>% select(Sample_ID,grep('s__',colnames(abundance_data)))

microbiome_vars = abundance_data %>% select(-Sample_ID) %>% colnames

# run log transform if not looking at deltas
if(!grepl('delta',abundance_data_file) & !grepl('diversity',abundance_data_file)){
	abundance_data = abundance_data %>% mutate_if(is.numeric,function(x) log(x+0.00001))
}

merged_data = inner_join(abundance_data, metadata, by='Sample_ID')

# compute associations of form microbe ~ age + treatment
regression_output_microbe_treatment = map(microbiome_vars, function(x) glm(data = merged_data, get(x) ~ age + rx) %>% tidy %>% filter(term!='(Intercept)',term!='age') %>% mutate(term = x)) %>% bind_rows  %>% mutate(bh = p.adjust(p.value))

saveRDS(regression_output_microbe_treatment,paste('metaphlan_associations/regression_output_microbe_treatment_',abundance_data_file,sep=''))