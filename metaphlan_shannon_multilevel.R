# metaphlan instead of bracken

library(tidyverse)
library(ggplot2)
library(umap)
library(cowplot)
library(vegan)
library(reshape2)
library(ggrepel)
library(broom)
library(ggpubr)

theme_set(theme_cowplot())

my_comparisons <- list( c("NOCHANGE-PLACEBO", "IMPROVED-PLACEBO"), c("NOCHANGE-TREATMENT", "RESPONDER"))

metadata = readRDS('pds08_metadata.rds') #%>% filter(b_bm_weekly<=4.2)
#demarcate responder vs nonresponder
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(rx != 'Treatment',"PLACEBO",if_else(d_bm_weekly>=1,'RESPONDER',if_else(abs(d_bm_weekly)<1,'NOCHANGE-TREATMENT','NOCHANGE-TREATMENT')))))
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(RESP_STATUS != 'PLACEBO',as.character(RESP_STATUS),if_else(d_bm_weekly>=1,'IMPROVED-PLACEBO',if_else(abs(d_bm_weekly)<1,'NOCHANGE-PLACEBO','NOCHANGE-PLACEBO')))))
metadata = metadata %>% mutate(IMPROVED = if_else(RESP_STATUS == 'IMPROVED-PLACEBO' | RESP_STATUS == 'RESPONDER','IMPROVED',if_else(RESP_STATUS == 'NOCHANGE-PLACEBO' | RESP_STATUS == 'NOCHANGE-TREATMENT','NOCHANGE','NOCHANGE')))

metadata = metadata %>% mutate(RESP_V_NONRESP = if_else(RESP_STATUS == 'RESPONDER',1,if_else(RESP_STATUS == 'NOCHANGE-TREATMENT',0,-1))) 
metadata$RESP_V_NONRESP[metadata$RESP_V_NONRESP==-1] = NA


baseline_metaphlan = readRDS('metaphlan_baseline.rds') %>% mutate(Sample_ID = paste(Sample_ID,'_baseline',sep=''))
endpoint_metaphlan = readRDS('metaphlan_endpoint.rds') %>% mutate(Sample_ID = paste(Sample_ID,'_endpoint',sep=''))
delta_metaphlan = readRDS('metaphlan_delta.rds') %>% mutate(Sample_ID = paste(Sample_ID,'_delta',sep=''))

data = bind_rows(baseline_metaphlan,endpoint_metaphlan) %>% column_to_rownames('Sample_ID') 
data = data %>% t

toremove = colnames(data) %>% strsplit('_') %>% map_chr(1) %>% table %>% data.frame %>% filter(Freq!=2) %>% select(colnames(.)[1]) %>% unlist %>% unname
data = data %>% data.frame%>% select(!(grep(paste(toremove, collapse="|"), colnames(data))))

#species
data_sub = data[grepl('s__',rownames(data)),]
shannon = map(colnames(data_sub), function(x) diversity(data_sub[,x],index='shannon')) %>% unlist %>% unname

divs = data.frame(shannon)
rownames(divs) = colnames(data)
divs = divs %>% rownames_to_column('Sample_ID')
colnames(divs) = c('Sample_ID','species')
divs_s = divs %>% mutate(timepoint = strsplit(Sample_ID,'_')%>% map_chr(2),Sample_ID = strsplit(Sample_ID,'_')%>% map_chr(1)) 

#genus
data_sub = data[!grepl('s__',rownames(data)),]
data_sub = data_sub[grepl('g__',rownames(data_sub)),]
shannon = map(colnames(data_sub), function(x) diversity(data_sub[,x],index='shannon')) %>% unlist %>% unname

divs = data.frame(shannon)
rownames(divs) = colnames(data)
divs = divs %>% rownames_to_column('Sample_ID')
colnames(divs) = c('Sample_ID','genus')
divs_g = divs %>% mutate(timepoint = strsplit(Sample_ID,'_')%>% map_chr(2),Sample_ID = strsplit(Sample_ID,'_')%>% map_chr(1)) 

#family
data_sub = data[!grepl('g__',rownames(data)),]
data_sub = data_sub[grepl('f__',rownames(data_sub)),]
shannon = map(colnames(data_sub), function(x) diversity(data_sub[,x],index='shannon')) %>% unlist %>% unname

divs = data.frame(shannon)
rownames(divs) = colnames(data)
divs = divs %>% rownames_to_column('Sample_ID')
colnames(divs) = c('Sample_ID','family')
divs_f = divs %>% mutate(timepoint = strsplit(Sample_ID,'_')%>% map_chr(2),Sample_ID = strsplit(Sample_ID,'_')%>% map_chr(1)) 

#order
data_sub = data[!grepl('f__',rownames(data)),]
data_sub = data_sub[grepl('o__',rownames(data_sub)),]
shannon = map(colnames(data_sub), function(x) diversity(data_sub[,x],index='shannon')) %>% unlist %>% unname

divs = data.frame(shannon)
rownames(divs) = colnames(data)
divs = divs %>% rownames_to_column('Sample_ID')
colnames(divs) = c('Sample_ID','order')
divs_o = divs %>% mutate(timepoint = strsplit(Sample_ID,'_')%>% map_chr(2),Sample_ID = strsplit(Sample_ID,'_')%>% map_chr(1)) 

#class
data_sub = data[!grepl('o__',rownames(data)),]
data_sub = data_sub[grepl('c__',rownames(data_sub)),]
shannon = map(colnames(data_sub), function(x) diversity(data_sub[,x],index='shannon')) %>% unlist %>% unname

divs = data.frame(shannon)
rownames(divs) = colnames(data)
divs = divs %>% rownames_to_column('Sample_ID')
colnames(divs) = c('Sample_ID','class')
divs_c = divs %>% mutate(timepoint = strsplit(Sample_ID,'_')%>% map_chr(2),Sample_ID = strsplit(Sample_ID,'_')%>% map_chr(1)) 

#phylum
data_sub = data[!grepl('c__',rownames(data)),]
data_sub = data_sub[grepl('p__',rownames(data_sub)),]
shannon = map(colnames(data_sub), function(x) diversity(data_sub[,x],index='shannon')) %>% unlist %>% unname

divs = data.frame(shannon)
rownames(divs) = colnames(data)
divs = divs %>% rownames_to_column('Sample_ID')
colnames(divs) = c('Sample_ID','phylum')
divs_p = divs %>% mutate(timepoint = strsplit(Sample_ID,'_')%>% map_chr(2),Sample_ID = strsplit(Sample_ID,'_')%>% map_chr(1)) 

# merge all data into one frame
#metaphlan_s = data %>% t %>% data.frame %>% select(all_of(grep('s__',colnames(.))))%>% rownames_to_column('Sample_ID') 
metaphlan_s = data %>% t %>% data.frame %>% rownames_to_column('Sample_ID') 
metaphlan_s$timepoint = strsplit(metaphlan_s$Sample_ID,'_')%>% map_chr(2)
metaphlan_s$Sample_ID = strsplit(metaphlan_s$Sample_ID,'_')%>% map_chr(1)
metaphlan_merged = left_join(metaphlan_s,divs_s,by=c('Sample_ID','timepoint'))
metaphlan_merged = left_join(metaphlan_merged,divs_g,by=c('Sample_ID','timepoint'))
metaphlan_merged = left_join(metaphlan_merged,divs_f,by=c('Sample_ID','timepoint'))
metaphlan_merged = left_join(metaphlan_merged,divs_o,by=c('Sample_ID','timepoint'))
metaphlan_merged = left_join(metaphlan_merged,divs_c,by=c('Sample_ID','timepoint'))
metaphlan_merged = left_join(metaphlan_merged,divs_p,by=c('Sample_ID','timepoint'))
metaphlan_merged = metaphlan_merged %>% select(-all_of(grep('k__',colnames(metaphlan_merged))))
metadata_sub = metadata %>% select(Sample_ID,RESP_STATUS)
abdata_div_melted = left_join(metaphlan_merged,metadata_sub,by='Sample_ID') %>% filter(!is.na(RESP_STATUS)) %>% melt

abdata_div_melted$RESP_STATUS = factor(abdata_div_melted$RESP_STATUS,levels = c("NOCHANGE-PLACEBO","IMPROVED-PLACEBO","NOCHANGE-TREATMENT","RESPONDER"))

ggplot(data = abdata_div_melted, aes(x = as.factor(timepoint), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +stat_compare_means(paired = T,method='wilcox.test')+ geom_point(aes(fill=as.factor(timepoint),group=Sample_ID), position = position_dodge(0.2)) + geom_line(aes(group=Sample_ID),color='grey',position = position_dodge(0.2)) + facet_grid(cols=vars(RESP_STATUS),rows=vars(variable),scales='free')
ggsave('richness_div_analysis/rnr_div_shannon_allclades.pdf',width=15,height=20)

ggplot(data = abdata_div_melted %>% filter(timepoint == 'baseline'), aes(x = as.factor(RESP_STATUS), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_point(aes(group=Sample_ID), position = position_dodge(0.2)) + facet_grid(rows=vars(variable),scales='free')+ stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')
ggsave('richness_div_analysis/rnr_div_shannon_base_allclades.pdf',width=15,height=20)

ggplot(data = abdata_div_melted %>% filter(timepoint == 'endpoint'), aes(x = as.factor(RESP_STATUS), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_point(aes(group=Sample_ID), position = position_dodge(0.2)) + facet_grid(rows=vars(variable),scales='free')+ stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')
ggsave('richness_div_analysis/rnr_div_shannon_end_allclades.pdf',width=15,height=20)

compare_means(data = abdata_div_melted %>% filter(timepoint == 'baseline') %>% filter(variable=='species') ,value ~ RESP_STATUS,method='wilcox.test')
compare_means(data = abdata_div_melted %>% filter(timepoint == 'baseline') %>% filter(variable=='genus') ,value ~ RESP_STATUS,method='wilcox.test')
compare_means(data = abdata_div_melted %>% filter(timepoint == 'baseline') %>% filter(variable=='family') ,value ~ RESP_STATUS,method='wilcox.test')
compare_means(data = abdata_div_melted %>% filter(timepoint == 'baseline') %>% filter(variable=='order') ,value ~ RESP_STATUS,method='wilcox.test')
compare_means(data = abdata_div_melted %>% filter(timepoint == 'baseline') %>% filter(variable=='class') ,value ~ RESP_STATUS,method='wilcox.test')
compare_means(data = abdata_div_melted %>% filter(timepoint == 'baseline') %>% filter(variable=='phylum') ,value ~ RESP_STATUS,method='wilcox.test')

### adjust for age to confirm baseline difference

abdata_div_melted_rnr = abdata_div_melted %>% filter(RESP_STATUS == 'RESPONDER' | RESP_STATUS == 'NOCHANGE-TREATMENT') %>% mutate(RESP_STATUS = as.numeric(as.factor(as.character(RESP_STATUS)))-1)

abdata_div_melted_rnr = left_join(abdata_div_melted_rnr,metadata %>% select(Sample_ID,age,b_bm_weekly))

s = glm(data = abdata_div_melted_rnr %>% filter(variable == 'species',timepoint == 'baseline'), RESP_STATUS ~ age + b_bm_weekly+ value,family = 'binomial') %>% tidy %>% mutate(level = 'species')
g = glm(data = abdata_div_melted_rnr %>% filter(variable == 'genus',timepoint == 'baseline'), RESP_STATUS ~ age +b_bm_weekly + value,family = 'binomial') %>% tidy%>% mutate(level = 'genus')
f = glm(data = abdata_div_melted_rnr %>% filter(variable == 'family',timepoint == 'baseline'), RESP_STATUS ~ age +b_bm_weekly + value,family = 'binomial') %>% tidy%>% mutate(level = 'family')
o = glm(data = abdata_div_melted_rnr %>% filter(variable == 'order',timepoint == 'baseline'), RESP_STATUS ~ age +b_bm_weekly+ value,family = 'binomial') %>% tidy%>% mutate(level = 'order')
c = glm(data = abdata_div_melted_rnr %>% filter(variable == 'class',timepoint == 'baseline'), RESP_STATUS ~ age +b_bm_weekly+ value,family = 'binomial') %>% tidy%>% mutate(level = 'class')
p = glm(data = abdata_div_melted_rnr %>% filter(variable == 'phylum',timepoint == 'baseline'), RESP_STATUS ~ age +b_bm_weekly+ value,family = 'binomial') %>% tidy%>% mutate(level = 'phylum')

glm_output = bind_rows(p,c,o,f,g,s)
write.csv(glm_output,'~/Desktop/google_drive/My Drive/pds08/richness_div_analysis/glm_output_metaphlan_shannon.csv')




