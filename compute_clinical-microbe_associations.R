#!/usr/bin/Rscript --vanilla

# compute associations with phenotypes

suppressMessages(library(tidyverse))
suppressMessages(library(broom))

#setwd('/Volumes/GoogleDrive-109433674545306273960/My\ Drive/pds08')
#setwd('~/Desktop/gdrive/pds08')

args <- commandArgs(trailingOnly = TRUE)

metadata_file = args[[1]]
abundance_data_file = args[[2]]
clinical_var = args[[3]]
datatype = args[[4]]

#metadata_file='pds08_metadata.rds'
#abundance_data_file = 'metaphlan_endpoint_diversity.rds'
#clinical_var = 'd_bm3'

metadata = readRDS(metadata_file)  %>% filter(b_bm_weekly<=4.2) %>% mutate(b_bm_weekly = log(b_bm_weekly)) 
abundance_data = readRDS(abundance_data_file)

# remove rows that are identical to each other and therefore encode no additional information + add hypotheses (big problem for metaphlan)
# this should take the HIGHER order
abundance_data = abundance_data %>% column_to_rownames('Sample_ID')
abundance_data = abundance_data[!duplicated(as.list(abundance_data))]
abundance_data = abundance_data %>% rownames_to_column('Sample_ID')

#abundance_data = abundance_data %>% select(Sample_ID,grep('s__',colnames(abundance_data)))

microbiome_vars = abundance_data %>% select(-Sample_ID) %>% colnames

# run log transform if not looking at deltas or diversity
if(!grepl('delta',abundance_data_file) & !grepl('diversity',abundance_data_file)){
	abundance_data = abundance_data %>% mutate_if(is.numeric,function(x) log(x+0.00001))
}

merged_data = inner_join(abundance_data, metadata, by='Sample_ID')

# compute associations of form outcome ~ age + microbe

family = 'gaussian'
if(length(unique(merged_data[,clinical_var])) - sum(is.na(merged_data[,clinical_var])) == 2L){
	family = 'binomial'
	merged_data[,clinical_var] = as.factor(merged_data[,clinical_var])
}

regression_output_outcome_microbe = map(microbiome_vars, function(x) glm(data = merged_data, family=family, get(clinical_var) ~ age + get(x)) %>% tidy %>% filter(term!='(Intercept)',term!='age') %>% mutate(dependent_var = clinical_var,term = x)) %>% bind_rows %>% mutate(bh = p.adjust(p.value))

saveRDS(regression_output_outcome_microbe,paste(datatype,'_associations/regression_output_outcome_microbe_',clinical_var,'_',abundance_data_file,sep=''))

