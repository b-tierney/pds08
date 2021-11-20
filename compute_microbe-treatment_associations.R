#!/usr/bin/Rscript --vanilla

# compute associations with phenotypes

suppressMessages(library(tidyverse))
suppressMessages(library(broom))

#setwd('/Volumes/GoogleDrive-109433674545306273960/My\ Drive/pds08')
#setwd('~/Desktop/gdrive/pds08')

args <- commandArgs(trailingOnly = TRUE)

metadata_file = args[[1]]
abundance_data_file = args[[2]]
datatype = args[[3]]

metadata = readRDS(metadata_file)
abundance_data = readRDS(abundance_data_file)

metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(rx != 'Treatment',"PLACEBO",if_else(d_bm_weekly>=1,'RESPONDER',if_else(abs(d_bm_weekly)<1,'NOCHANGE-TREATMENT','NOCHANGE-TREATMENT')))))
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(RESP_STATUS != 'PLACEBO',as.character(RESP_STATUS),if_else(d_bm_weekly>=1,'IMPROVED-PLACEBO',if_else(abs(d_bm_weekly)<1,'NOCHANGE-PLACEBO','NOCHANGE-PLACEBO')))))
metadata = metadata %>% mutate(IMPROVED = if_else(RESP_STATUS == 'IMPROVED-PLACEBO' | RESP_STATUS == 'RESPONDER','IMPROVED',if_else(RESP_STATUS == 'NOCHANGE-PLACEBO' | RESP_STATUS == 'NOCHANGE-TREATMENT','NOCHANGE','NOCHANGE')))

metadata = metadata %>% mutate(RESP_V_NONRESP = if_else(RESP_STATUS == 'RESPONDER',1,if_else(RESP_STATUS == 'NOCHANGE-TREATMENT',0,-1))) 
metadata$RESP_V_NONRESP[metadata$RESP_V_NONRESP==-1] = NA

# remove rows that are identical to each other and therefore encode no additional information + add hypotheses (big problem for metaphlan)
# this should take the HIGHER order
abundance_data = abundance_data %>% column_to_rownames('Sample_ID')
abundance_data = abundance_data[!duplicated(as.list(abundance_data))]
abundance_data = abundance_data %>% rownames_to_column('Sample_ID')

#abundance_data = abundance_data %>% select(Sample_ID,grep('s__',colnames(abundance_data)))

microbiome_vars = abundance_data %>% select(-Sample_ID) %>% colnames

# run log transform if not looking at deltas
if(!grepl('delta',abundance_data_file) & !grepl('diversity',abundance_data_file)){
	abundance_data = abundance_data %>% mutate_if(is.numeric,function(x) log(x+0.00001))
}

merged_data = inner_join(abundance_data, metadata, by='Sample_ID')
merged_data = merged_data %>% filter(!is.na(RESP_V_NONRESP))

# compute associations of form microbe ~ age + treatment
regression_output_microbe_treatment = map(microbiome_vars, function(x) glm(data = merged_data, family = 'binomial',RESP_V_NONRESP ~ age + b_bm_weekly + get(x)) %>% tidy %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly') %>% mutate(term = x)) %>% bind_rows  %>% mutate(bh = p.adjust(p.value))

saveRDS(regression_output_microbe_treatment,paste(datatype,'_associations/regression_output_microbe_treatment_',abundance_data_file,sep=''))









