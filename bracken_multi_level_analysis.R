
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


setwd('~/Desktop/google_drive/My Drive/SEED_HEALTH/PDS08 (1)/')

build_bracken_abmat <- function(brackenfile,level){
  bracken = read.csv(brackenfile,sep='\t',header=T)
  bracken_num = bracken %>% select(name,all_of(colnames(bracken)[!grepl('frac',colnames(bracken))])) %>% column_to_rownames('name') %>% select(-taxonomy_id,-taxonomy_lvl)
  orgs_to_keep = rowMeans(bracken_num) %>% data.frame %>% rownames_to_column('organism') %>% select(organism) %>% unlist %>% unname
  
  bracken = bracken %>% select(name,taxonomy_lvl,all_of(colnames(bracken)[grepl('frac',colnames(bracken))]))  %>% filter(name %in% orgs_to_keep)
  colnames(bracken) = gsub('\\.bracken_frac','',colnames(bracken))
  
  # get samples to keep 
  basekeep = colnames(bracken)[grepl('_01',colnames(bracken))] 
  basekeep = gsub('_01','',basekeep)
  endkeep = colnames(bracken)[grepl('_02',colnames(bracken))] 
  endkeep = gsub('_02','',endkeep)
  
  samples_to_keep = intersect(basekeep,endkeep)
  
  colnames(bracken) = gsub('_01','_baseline',colnames(bracken))
  colnames(bracken) = gsub('_02','_endpoint',colnames(bracken))
  colnames(bracken) = gsub('_family','',colnames(bracken))
  
  bracken_f = bracken %>% select(-taxonomy_lvl)
  my_comparisons <- list( c("NOCHANGE-PLACEBO", "IMPROVED-PLACEBO"), c("NOCHANGE-TREATMENT", "RESPONDER"))
  
  metadata = readRDS('pds08_metadata.rds') #%>% filter(b_bm_weekly<=4.2)
  #demarcate responder vs nonresponder
  metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(rx != 'Treatment',"PLACEBO",if_else(d_bm_weekly>=1,'RESPONDER',if_else(abs(d_bm_weekly)<1,'NOCHANGE-TREATMENT','NOCHANGE-TREATMENT')))))
  metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(RESP_STATUS != 'PLACEBO',as.character(RESP_STATUS),if_else(d_bm_weekly>=1,'IMPROVED-PLACEBO',if_else(abs(d_bm_weekly)<1,'NOCHANGE-PLACEBO','NOCHANGE-PLACEBO')))))
  metadata = metadata %>% mutate(IMPROVED = if_else(RESP_STATUS == 'IMPROVED-PLACEBO' | RESP_STATUS == 'RESPONDER','IMPROVED',if_else(RESP_STATUS == 'NOCHANGE-PLACEBO' | RESP_STATUS == 'NOCHANGE-TREATMENT','NOCHANGE','NOCHANGE')))
  
  metadata = metadata %>% mutate(RESP_V_NONRESP = if_else(RESP_STATUS == 'RESPONDER',1,if_else(RESP_STATUS == 'NOCHANGE-TREATMENT',0,-1))) 
  metadata$RESP_V_NONRESP[metadata$RESP_V_NONRESP==-1] = NA
  
  #compute diversity and richness
  data = bracken_f %>% column_to_rownames('name')
  
  toremove = colnames(data) %>% strsplit('_') %>% map_chr(1) %>% table %>% data.frame %>% filter(Freq!=2) %>% select(colnames(.)[1]) %>% unlist %>% unname
  data = data %>% data.frame%>% select(!(grep(paste(toremove, collapse="|"), colnames(data))))
  
  #data = data[rowMeans(data)>.00001,] 
  shannon = map(colnames(data), function(x) diversity(data[,x],index='shannon')) %>% unlist %>% unname
  simpson = map(colnames(data), function(x) diversity(data[,x],index='simpson')) %>% unlist %>% unname
  richness = map(colnames(data), function(x) length(which(data[,x]!=0))) %>% unlist %>% unname
  
  divs = data.frame(shannon,simpson,richness)
  rownames(divs) = colnames(data)
  divs = divs %>% rownames_to_column('Sample_ID')
  divs = divs %>% mutate(timepoint = strsplit(Sample_ID,'_')%>% map_chr(2),Sample_ID = strsplit(Sample_ID,'_')%>% map_chr(1)) 
  
  # merge all data into one frame
  bracken_f = bracken_f %>% column_to_rownames('name') %>% t %>% data.frame %>% rownames_to_column('Sample_ID') 
  bracken_f$timepoint = strsplit(bracken_f$Sample_ID,'_')%>% map_chr(2)
  bracken_f$Sample_ID = strsplit(bracken_f$Sample_ID,'_')%>% map_chr(1)
  
  bracken_merged = left_join(bracken_f,divs,by=c('Sample_ID','timepoint'))
  abdata = left_join(bracken_merged,metadata,by='Sample_ID') %>% filter(!is.na(RESP_STATUS))

  abdata_div = abdata %>%  select(Sample_ID,shannon,simpson,richness,RESP_STATUS,rx,timepoint,IMPROVED)
  abdata_div_melted = melt(abdata_div)
  
  abdata_div_melted$RESP_STATUS = factor(abdata_div_melted$RESP_STATUS,levels = c("NOCHANGE-PLACEBO","IMPROVED-PLACEBO","NOCHANGE-TREATMENT","RESPONDER"))
  
  ggplot(data = abdata_div_melted, aes(x = as.factor(timepoint), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +stat_compare_means(paired = T,method='wilcox.test')+ geom_point(aes(fill=as.factor(timepoint),group=Sample_ID), position = position_dodge(0.2)) + geom_line(aes(group=Sample_ID),color='grey',position = position_dodge(0.2)) + facet_grid(cols=vars(RESP_STATUS),rows=vars(variable),scales='free')
  ggsave(paste('richness_div_analysis_bracken/rnr_div_',level,'.pdf',sep=''),width=15,height=15)
  
  ggplot(data = abdata_div_melted %>% filter(timepoint == 'baseline'), aes(x = as.factor(RESP_STATUS), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_point(aes(group=Sample_ID), position = position_dodge(0.2)) + facet_grid(rows=vars(variable),scales='free')+ stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')
  ggsave(paste('richness_div_analysis_bracken/rnr_div_base_',level,'.pdf',sep=''),width=15,height=15)
  
  ggplot(data = abdata_div_melted %>% filter(timepoint == 'endpoint'), aes(x = as.factor(RESP_STATUS), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_point(aes(group=Sample_ID), position = position_dodge(0.2)) + facet_grid(rows=vars(variable),scales='free')+ stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')
  ggsave(paste('richness_div_analysis_bracken/rnr_div_end_',level,'.pdf',sep=''),width=15,height=15)
  
  toreturn = ggplot(data = abdata_div_melted %>% filter(timepoint == 'baseline',variable == 'richness'), aes(x = as.factor(RESP_STATUS), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_point(aes(group=Sample_ID), position = position_dodge(0.2)) + ggtitle(toupper(level)) + stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')
  toreturn2 = abdata_div_melted %>% mutate(level=level)
  
  abdata_div_melted_rnr = abdata_div_melted %>% filter(RESP_STATUS == 'RESPONDER' | RESP_STATUS == 'NOCHANGE-TREATMENT') %>% mutate(RESP_STATUS = as.numeric(as.factor(as.character(RESP_STATUS)))-1)
  
  abdata_div_melted_rnr = left_join(abdata_div_melted_rnr,metadata %>% select(Sample_ID,age,b_bm_weekly))
  
  glmout = glm(data = abdata_div_melted_rnr %>% filter(timepoint == 'baseline',variable == 'richness'), RESP_STATUS ~ age +b_bm_weekly+ value,family = 'binomial') %>% tidy %>% mutate(level = level)
  return(list(toreturn,glmout,toreturn2))
}

species = build_bracken_abmat("bracken_output/merged_bracken_output.tsv",'species')
ggsave(plot = species[[1]],'richness_div_analysis_bracken/species_singleplot.pdf',width=8,height=8)
genus = build_bracken_abmat("bracken_output/merged_bracken_output_genus.tsv",'genus')
ggsave(plot = genus[[1]],'richness_div_analysis_bracken/genus_singleplot.pdf',width=8,height=8)
family = build_bracken_abmat("bracken_output/merged_bracken_output_families.tsv",'family')
ggsave(plot = family[[1]],'richness_div_analysis_bracken/family_singleplot.pdf',width=8,height=8)
order = build_bracken_abmat("bracken_output/merged_bracken_output_orders.tsv",'order')
ggsave(plot = order[[1]],'richness_div_analysis_bracken/order_singleplot.pdf',width=8,height=8)
class = build_bracken_abmat("bracken_output/merged_bracken_output_class.tsv",'class')
ggsave(plot = class[[1]],'richness_div_analysis_bracken/class_singleplot.pdf',width=8,height=8)
phylum = build_bracken_abmat("bracken_output/merged_bracken_output_phyla.tsv",'phylum')
ggsave(plot = phylum[[1]],'richness_div_analysis_bracken/phylum_singleplot.pdf',width=8,height=8)

glm_output = bind_rows(phylum[[2]],class[[2]],order[[2]],family[[2]],genus[[2]],species[[2]])
write.csv(glm_output,'~/Desktop/google_drive/My Drive/pds08/richness_div_analysis_bracken//glm_output_bracken.csv')

rawdata = bind_rows(phylum[[3]],class[[3]],order[[3]],family[[3]],genus[[3]],species[[3]])



