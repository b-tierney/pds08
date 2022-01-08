#!/usr/bin/Rscript --vanilla

library(tidyverse)

#setwd('~/Desktop/google_drive/My Drive/pds08')

### metadata cleaning

# load data
metadata = read.csv('meta.csv')
sequencing_data = read.csv('sequencing_metadata.csv') %>% select(Seed_ID,Curebase_ID,Cryo_Vial_Serial_Number,Randomization_Assignment) %>% rename(Sample_ID = Cryo_Vial_Serial_Number)
#sequencing_data2 = read.csv('sequencing_metadata_2.csv') %>% select(Sample_ID,Timepoint) %>% mutate(Sample_ID = strsplit(Sample_ID,'_') %>% map_chr(1))

merged_metadata = left_join(metadata,sequencing_data,by=c('id'='Curebase_ID')) %>% filter(!is.na(Sample_ID))

merged_metadata = merged_metadata %>% mutate(rx = if_else(rx == 'Active','Treatment',rx))

#demarcate responder vs nonresponder
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(rx != 'Treatment',"PLACEBO",if_else(d_bm_weekly>=1,'RESPONDER',if_else(abs(d_bm_weekly)<1,'NOCHANGE-TREATMENT','NOCHANGE-TREATMENT')))))
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(RESP_STATUS != 'PLACEBO',as.character(RESP_STATUS),if_else(d_bm_weekly>=1,'IMPROVED-PLACEBO',if_else(abs(d_bm_weekly)<1,'NOCHANGE-PLACEBO','NOCHANGE-PLACEBO')))))
metadata = metadata %>% mutate(IMPROVED = if_else(RESP_STATUS == 'IMPROVED-PLACEBO' | RESP_STATUS == 'RESPONDER','IMPROVED',if_else(RESP_STATUS == 'NOCHANGE-PLACEBO' | RESP_STATUS == 'NOCHANGE-TREATMENT','NOCHANGE','NOCHANGE')))

metadata = metadata %>% mutate(RESP_V_NONRESP = if_else(RESP_STATUS == 'RESPONDER',1,if_else(RESP_STATUS == 'NOCHANGE-TREATMENT',0,-1))) 
metadata$RESP_V_NONRESP[metadata$RESP_V_NONRESP==-1] = NA
saveRDS(merged_metadata,'pds08_metadata.rds')
### abundance cleaning

# metaphlan
metaphlan = read.table('pds08_metaphlan_output.txt',header=T) %>% select(-NCBI_tax_id)
metaphlan_baseline = metaphlan %>% select(clade_name,grep('_01',colnames(metaphlan))) %>% column_to_rownames('clade_name')
colnames(metaphlan_baseline) = map(colnames(metaphlan_baseline), function(x) gsub('_01_metaphlan3_out','',x))
metaphlan_endpoint = metaphlan %>% select(clade_name,grep('_02',colnames(metaphlan))) %>% column_to_rownames('clade_name')
colnames(metaphlan_endpoint) = map(colnames(metaphlan_endpoint), function(x) gsub('_02_metaphlan3_out','',x))

matched_columns = intersect(colnames(metaphlan_baseline),colnames(metaphlan_endpoint))

metaphlan_baseline_sub = metaphlan_baseline %>% select(all_of(matched_columns))
metaphlan_endpoint_sub = metaphlan_endpoint %>% select(all_of(matched_columns))

metaphlan_deltas = metaphlan_endpoint_sub - metaphlan_baseline_sub

metaphlan_baseline = t(metaphlan_baseline) %>% as.data.frame%>% rownames_to_column('Sample_ID')
metaphlan_endpoint = t(metaphlan_endpoint)%>% as.data.frame %>% rownames_to_column('Sample_ID')
metaphlan_deltas = t(metaphlan_deltas) %>% as.data.frame %>% rownames_to_column('Sample_ID')

saveRDS(metaphlan_baseline,'metaphlan_baseline.rds')
saveRDS(metaphlan_endpoint,'metaphlan_endpoint.rds')
saveRDS(metaphlan_deltas,'metaphlan_delta.rds')

# pathways

pathways = read.csv('pds08_pathways_output.txt',sep='\t')

colnames(pathways)[1]='pathway'

pathways = pathways[pathways$pathway!='UNMAPPED',]
pathways = pathways[!grepl('UNINTEGRATED',pathways$pathway),] 
pathways = pathways[!grepl('\\|',pathways$pathway),] 

pathways_baseline = pathways %>% select(pathway,grep('_01',colnames(pathways))) %>% as_tibble %>% column_to_rownames('pathway')
colnames(pathways_baseline) = map(colnames(pathways_baseline), function(x) gsub('_01_all_reads_Abundance','',x))
pathways_endpoint = pathways %>% select(pathway,grep('_02',colnames(pathways)))  %>% as_tibble %>% column_to_rownames('pathway')
colnames(pathways_endpoint) = map(colnames(pathways_endpoint), function(x) gsub('_02_all_reads_Abundance','',x))

matched_columns = intersect(colnames(pathways_baseline),colnames(pathways_endpoint))

pathways_baseline_sub = pathways_baseline %>% select(all_of(matched_columns))
pathways_endpoint_sub = pathways_endpoint %>% select(all_of(matched_columns))

pathways_deltas = pathways_endpoint_sub - pathways_baseline_sub

pathways_baseline = t(pathways_baseline) %>% as.data.frame%>% rownames_to_column('Sample_ID')
pathways_endpoint = t(pathways_endpoint)%>% as.data.frame %>% rownames_to_column('Sample_ID')
pathways_deltas = t(pathways_deltas) %>% as.data.frame %>% rownames_to_column('Sample_ID')

saveRDS(pathways_baseline,'pathways_baseline.rds')
saveRDS(pathways_endpoint,'pathways_endpoint.rds')
saveRDS(pathways_deltas,'pathways_delta.rds')

# gene families
genefamilies = read.csv('pds08_genefamilies_output.txt',sep='\t')

colnames(genefamilies)[1]='gene_family'

genefamilies = genefamilies[genefamilies$gene_family!='UNMAPPED',]
genefamilies = genefamilies[!grepl('unclassified',genefamilies$gene_family),] 
genefamilies = genefamilies[!grepl('\\|',genefamilies$gene_family),] 

genefamilies = genefamilies %>% as_tibble %>% column_to_rownames('gene_family')

genefamilies = genefamilies[rowMeans(genefamilies == 0) <= 0.9,] %>% rownames_to_column('gene_family')

genefamilies_baseline = genefamilies %>% select(gene_family,grep('_01',colnames(genefamilies))) %>% as_tibble %>% column_to_rownames('gene_family')
colnames(genefamilies_baseline) = map(colnames(genefamilies_baseline), function(x) gsub('_01_all_reads_Abundance.RPKs','',x))
genefamilies_endpoint = genefamilies %>% select(gene_family,grep('_02',colnames(genefamilies)))  %>% as_tibble %>% column_to_rownames('gene_family')
colnames(genefamilies_endpoint) = map(colnames(genefamilies_endpoint), function(x) gsub('_02_all_reads_Abundance.RPKs','',x))

matched_columns = intersect(colnames(genefamilies_baseline),colnames(genefamilies_endpoint))

genefamilies_baseline_sub = genefamilies_baseline %>% select(all_of(matched_columns))
genefamilies_endpoint_sub = genefamilies_endpoint %>% select(all_of(matched_columns))

genefamilies_deltas = genefamilies_endpoint_sub - genefamilies_baseline_sub

genefamilies_baseline = t(genefamilies_baseline) %>% as.data.frame%>% rownames_to_column('Sample_ID')
genefamilies_endpoint = t(genefamilies_endpoint)%>% as.data.frame %>% rownames_to_column('Sample_ID')
genefamilies_deltas = t(genefamilies_deltas) %>% as.data.frame %>% rownames_to_column('Sample_ID')

saveRDS(genefamilies_baseline,'genefamilies_baseline.rds')
saveRDS(genefamilies_endpoint,'genefamilies_endpoint.rds')
saveRDS(genefamilies_deltas,'genefamilies_delta.rds')










