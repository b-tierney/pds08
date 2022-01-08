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

setwd('~/Desktop/google_drive/My Drive/pds08/')

metadata = readRDS('pds08_metadata.rds') 
#demarcate responder vs nonresponder
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(rx != 'Treatment',"PLACEBO",if_else(d_bm_weekly>=1,'RESPONDER',if_else(abs(d_bm_weekly)<1,'NOCHANGE-TREATMENT','NOCHANGE-TREATMENT')))))
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(RESP_STATUS != 'PLACEBO',as.character(RESP_STATUS),if_else(d_bm_weekly>=1,'IMPROVED-PLACEBO',if_else(abs(d_bm_weekly)<1,'NOCHANGE-PLACEBO','NOCHANGE-PLACEBO')))))
metadata = metadata %>% mutate(IMPROVED = if_else(RESP_STATUS == 'IMPROVED-PLACEBO' | RESP_STATUS == 'RESPONDER','IMPROVED',if_else(RESP_STATUS == 'NOCHANGE-PLACEBO' | RESP_STATUS == 'NOCHANGE-TREATMENT','NOCHANGE','NOCHANGE')))

metadata = metadata %>% mutate(RESP_V_NONRESP = if_else(RESP_STATUS == 'RESPONDER',1,if_else(RESP_STATUS == 'NOCHANGE-TREATMENT',0,-1))) 
metadata$RESP_V_NONRESP[metadata$RESP_V_NONRESP==-1] = NA


metadata$rx = as.numeric(as.factor(metadata$rx))-1

baseline_metaphlan = readRDS('metaphlan_baseline.rds') %>% mutate(Sample_ID = paste(Sample_ID,'_baseline',sep=''))
endpoint_metaphlan = readRDS('metaphlan_endpoint.rds') %>% mutate(Sample_ID = paste(Sample_ID,'_endpoint',sep=''))

data = bind_rows(baseline_metaphlan,endpoint_metaphlan) %>% column_to_rownames('Sample_ID') 
data = data %>% t

toremove = colnames(data) %>% strsplit('_') %>% map_chr(1) %>% table %>% data.frame %>% filter(Freq!=2) %>% select(colnames(.)[1]) %>% unlist %>% unname
data = data %>% data.frame%>% select(!(grep(paste(toremove, collapse="|"), colnames(data))))

# merge all data into one frame
#metaphlan_s = data %>% t %>% data.frame %>% select(all_of(grep('s__',colnames(.))))%>% rownames_to_column('Sample_ID') 
metaphlan_s = data %>% t %>% data.frame %>% rownames_to_column('Sample_ID') 
metaphlan_s$timepoint = strsplit(metaphlan_s$Sample_ID,'_')%>% map_chr(2)
metaphlan_s$Sample_ID = strsplit(metaphlan_s$Sample_ID,'_')%>% map_chr(1)
abdata = left_join(metaphlan_s,metadata,by='Sample_ID') %>% filter(!is.na(RESP_STATUS))

######### ALL CLINICAL ASSOCIATIONS, ONLY USE METADATA FILTERED < BM5 AND PREVALENCE > 3

base_prevs = abdata %>% filter(timepoint == 'baseline',b_bm_weekly<5) %>% select(all_of(grep('k__',colnames(abdata %>% filter(timepoint == 'baseline')))))
base_prevs[base_prevs>0]=1
base_prevs = colSums(base_prevs) %>% data.frame %>% rownames_to_column('organism')
colnames(base_prevs) = c('organism','prev')
base_regcolumns = base_prevs %>% filter(prev>=10) %>% select(organism) %>% unlist %>% unname

end_prevs = abdata %>% filter(timepoint == 'endpoint',b_bm_weekly<5) %>% select(all_of(grep('k__',colnames(abdata %>% filter(timepoint == 'endpoint')))))
end_prevs[end_prevs>0]=1
end_prevs = colSums(end_prevs) %>% data.frame %>% rownames_to_column('organism')
colnames(end_prevs) = c('organism','prev')
end_regcolumns = end_prevs %>% filter(prev>=10) %>% select(organism) %>% unlist %>% unname

### TREATMENT ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , family = 'binomial', rx ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "treatment__endpoint-abundances")
}
output_merged_RX_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### D_BM2 ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , family = 'binomial', d_bm2 ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var') %>% mutate(regression = "dbm2__endpoint-abundances")
}
output_merged_DBM2_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### D_BM3 ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , family = 'binomial', d_bm3 ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "dbm3__endpoint-abundances")
}
output_merged_DBM3_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### B_BM_WEEKLY ~ BASELINE METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'baseline')

output = list()
for(r in base_regcolumns){
  regout = glm(data = regdata , b_bm_weekly ~ age + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "baseline-weekly-bms__baseline-abundances")
}
output_merged_B_BM_WEEKLY_BASELINE = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### BM_WEEKLY ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , bm_weekly ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "endpoint-weekly-bms__endpoint-abundances")
}
output_merged_BM_WEEKLY_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### B_BLOAT ~ BASELINE METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'baseline')

output = list()
for(r in base_regcolumns){
  regout = glm(data = regdata , b_bloat ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "baseline-bloating__baseline-abundances")
}
output_merged_B_BLOAT = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### D_BLOAT ~ BASELINE METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'baseline')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , d_bloat ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')
}
output_merged_D_BLOAT_BASE = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) %>% mutate(regression = "delta-bloating__baseline-abundances")


### D_BLOAT ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , d_bloat ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "delta-bloating__endpoint-abundances")
}
output_merged_D_BLOAT_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### BLOAT ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , bloat ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "endpoint-bloating__endpoint-abundances")
}
output_merged_BLOAT_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### B_PAIN ~ BASELINE METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'baseline')

output = list()
for(r in base_regcolumns){
  regout = glm(data = regdata , b_pain ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "baseline-pain__baseline-abundances")
}
output_merged_B_PAIN = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### D_PAIN ~ BASELINE METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'baseline')

output = list()
for(r in base_regcolumns){
  regout = glm(data = regdata , d_pain ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "delta-pain__baseline-abundances")
}
output_merged_D_PAIN_BASE = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 


### D_PAIN ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , d_pain ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var') %>% mutate(regression = "delta-pain__endpoint-abundances")
}
output_merged_D_PAIN_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### PAIN ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(b_bm_weekly<5) %>% filter(timepoint == 'endpoint')

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , pain ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var') %>% mutate(regression = "endpoint-pain__endpoint-abundances")
}
output_merged_PAIN_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

clinical_taxonomic_regressions = bind_rows(output_merged_RX_END, output_merged_D_BLOAT_END,output_merged_DBM2_END, output_merged_DBM3_END, output_merged_BM_WEEKLY_END, output_merged_B_BLOAT, output_merged_D_BLOAT_BASE, output_merged_D_BLOAT_BASE, output_merged_BLOAT_END, output_merged_B_PAIN, output_merged_D_PAIN_BASE, output_merged_D_PAIN_END,output_merged_PAIN_END)

regressions = unique(clinical_taxonomic_regressions$regression)

for(r in regressions){
  temp  = clinical_taxonomic_regressions %>% filter(regression == r) %>% mutate(FDR_SIGNIFICANT = if_else(adj<0.05,1,0))
  ggplot(data = temp,aes(x=estimate,color=factor(FDR_SIGNIFICANT),y=-log(p.value,10))) + geom_point(aes(alpha=.5)) + ylab('p-value') + geom_label_repel(data = temp %>% arrange(p.value) %>% filter(p.value<0.05) %>% arrange(p.value) %>% head(3),aes(label = microbial_var),box.padding   = 0.1, point.padding = 0.1,max.overlaps=100,size=2,segment.color = 'grey50')+geom_hline(yintercept = -log(0.05,5)) + xlim(-2.5,2.5) + ylim(0,6) + ggtitle(toupper(r))+ theme(legend.position = "none")
  ggsave(paste("./regression_output/",r,'.pdf',sep=''),width=6,height=6)
}

write.csv(clinical_taxonomic_regressions,'./regression_output/taxonomic_clinical_regressions.csv')
######### RESPONSE TO TREATMENT ASSOCIATIONS

### RESP_V_NONRESP ~ BASELINE METAPHLAN
regdata = abdata %>% filter(timepoint == 'baseline') %>% filter(!is.na(RESP_V_NONRESP))

output = list()
for(r in base_regcolumns){
  regout = glm(data = regdata , family = 'binomial', RESP_V_NONRESP ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "resp-v-nonresp-fullcohort__baseline-abundances")
}
output_merged_RESPVNONRESP_BASE = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

### RESP_V_NONRESP ~ ENDPOINT METAPHLAN
regdata = abdata %>% filter(timepoint == 'endpoint') %>% filter(!is.na(RESP_V_NONRESP))

output = list()
for(r in end_regcolumns){
  regout = glm(data = regdata , family = 'binomial', RESP_V_NONRESP ~ age + b_bm_weekly + log(regdata[,r]+0.00001)) %>% tidy %>% mutate(microbial_var = r) %>% filter(term!='(Intercept)',term!='age',term!='b_bm_weekly')
  prevs = regdata %>% select(r) %>% mutate(prev = if_else(get(r)!=0,1,0))  %>% summarize(prev=sum(prev)) %>% mutate(microbial_var = r)
  output[[r]] =left_join(regout,prevs,by='microbial_var')%>% mutate(regression = "resp-v-nonresp-fullcohort__endpoint-abundances")
}
output_merged_RESPVNONRESP_END = bind_rows(output) %>% mutate(adj=p.adjust(p.value,method = 'BY')) %>% select(-term) 

rnr_taxonomic = bind_rows(output_merged_RESPVNONRESP_BASE,output_merged_RESPVNONRESP_END)

write.csv(rnr_taxonomic,'./regression_output/taxonomic_responder_vs_non_responder_analysis.csv')