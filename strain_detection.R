library(tidyverse)

setwd('~/Desktop/google_drive/My Drive/SEED_HEALTH/PDS08 (1)/metapan/')

mapping = read.csv('metapangenome_contig_prefix_mapping',header=F,sep=' ') %>% mutate(V1 = strsplit(V1,'\\.') %>% map_chr(1))
colnames(mapping) = c('species','anvio_id')

files = list.files()
files = files[grepl('GENE-COVs',files)]

output1 = list()
output2=list()
for(f in files){
  data = read.csv(f,sep='\t') %>% select(-key)
  totalgenes = nrow(data)
  #data_base = data %>% select(all_of(grep('02',colnames(data))))
  data[data!=0]=1
  fract = colSums(data!=0)/totalgenes
  fract = as.data.frame(fract)
  fract = fract %>% rownames_to_column('id') %>% mutate(timepoint = strsplit(id,'_') %>% map_chr(2)) %>% mutate(timepoint = if_else(timepoint == '01','baseline','endpoint')) %>% mutate(Sample_ID = id %>% strsplit('_') %>% map_chr(1)) %>% select(-id) %>% mutate(base = strsplit(f,'-') %>% map_chr(1))
  output2[[f]] = fract
  data = log(read.csv(f,sep='\t') %>% select(-key)+0.00001) %>% colMeans() %>% data.frame %>% rownames_to_column('id') %>% mutate(timepoint = strsplit(id,'_') %>% map_chr(2)) %>% mutate(timepoint = if_else(timepoint == '01','baseline','endpoint')) %>% mutate(Sample_ID = id %>% strsplit('_') %>% map_chr(1)) %>% select(-id)
  base = strsplit(f,'-') %>% map_chr(1)
  colnames(data)[1] = mapping %>% filter(anvio_id == base) %>% select(species) %>% unlist %>% unname
  output1[[f]] = data %>% melt
}

output = bind_rows(output1)

metadata = readRDS('../pds08_metadata.rds')
#demarcate responder vs nonresponder
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(rx != 'Treatment',"PLACEBO",if_else(d_bm_weekly>=2,'RESPONDER',if_else(abs(d_bm_weekly)<2,'NOCHANGE-TREATMENT','WORSE-TREATMENT')))))
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(RESP_STATUS != 'PLACEBO',as.character(RESP_STATUS),if_else(d_bm_weekly>=2,'IMPROVED-PLACEBO',if_else(abs(d_bm_weekly)<2,'NOCHANGE-PLACEBO','WORSE-PLACEBO')))))
metadata = metadata %>% mutate(IMPROVED = if_else(RESP_STATUS == 'IMPROVED-PLACEBO' | RESP_STATUS == 'RESPONDER','IMPROVED',if_else(RESP_STATUS == 'NOCHANGE-PLACEBO' | RESP_STATUS == 'NOCHANGE-TREATMENT','NOCHANGE','WORSE')))

merged = left_join(output,metadata %>% select(Sample_ID,RESP_STATUS,rx),by='Sample_ID') %>% filter(!is.na(RESP_STATUS))

tokeep = merged %>% select(timepoint,Sample_ID) %>% unique %>% select(-timepoint) %>% table %>% data.frame %>% filter(Freq==2) %>% select(".") %>% unlist %>% unname %>% as.character

merged = merged %>% filter(Sample_ID %in% tokeep)

ggplot(data = merged, aes(x = as.factor(timepoint), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_point(aes(fill=as.factor(timepoint),group=Sample_ID), position = position_dodge(0.2)) + geom_line(aes(group=Sample_ID),color='grey',position = position_dodge(0.2)) + facet_grid(rows=vars(rx),cols=vars(variable),scales='free')  +stat_compare_means(paired = T,method='wilcox.test')
ggsave('./strain_analysis_rx.pdf',widt=25,height=20)

ggplot(data = merged %>% filter(timepoint == 'endpoint'), aes(x = as.factor(rx), y = value)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ geom_point(aes(group=Sample_ID), position = position_dodge(0.2)) + facet_grid(cols=vars(variable),scales='free')+stat_compare_means(paired = F,method='wilcox.test')
ggsave('strain_analysis_end_rx.pdf',width=25,height=7)

# get number of individuals with a given probiotic strain 
output2 = bind_rows(output2)
merged = left_join(output2,metadata %>% select(Sample_ID,RESP_STATUS,rx),by='Sample_ID') %>% filter(!is.na(RESP_STATUS))
tokeep = merged %>% select(timepoint,Sample_ID) %>% unique %>% select(-timepoint) %>% table %>% data.frame %>% filter(Freq==2) %>% select(".") %>% unlist %>% unname %>% as.character
merged = merged %>% filter(Sample_ID %in% tokeep)

merged = left_join(merged,mapping,by=c('base'='anvio_id'))

