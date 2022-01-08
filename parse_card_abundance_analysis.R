#parse card gene output

library(tidyverse)
library(broom)

files = read.csv('genemappinglocs',header=F) %>% unlist %>% unname

# load and normalize card bwt output

abundance_data = list()
for(f in files){
  data = read.csv(f,sep='\t') %>% select(ARO.Term,All.Mapped.Reads,Reference.Length,AMR.Gene.Family,Drug.Class,Resistance.Mechanism)
  if(nrow(data)>0){
  data$Reference.Length = as.character(data$Reference.Length)
  data = data %>% mutate(Reference.Length = Reference.Length %>% strsplit(';') %>% map_chr(1) %>% as.numeric,divfrac = as.numeric(data$All.Mapped.Reads)/Reference.Length, relative_abundance = divfrac/sum(divfrac))
  data = data %>% select(ARO.Term,relative_abundance) %>% mutate(id = f %>% strsplit('/') %>% map_chr(2) %>% gsub('_rgibwt_output','',.), timepoint = id %>% strsplit('_') %>% map_chr(2),timepoint = if_else(timepoint == '01','baseline','endpoint'),Sample_ID = id %>% strsplit('_') %>% map_chr(1)) %>% select(-id)
  abundance_data[[f]] = data
  }
}

abundance_data = bind_rows(abundance_data)

# get delta data

abundance_data_b = abundance_data %>% filter(timepoint == 'baseline') %>% select(-timepoint)
abundance_data_e = abundance_data %>% filter(timepoint == 'endpoint')%>% select(-timepoint)
delta_abundance = left_join(abundance_data_b,abundance_data_e,by = c('ARO.Term','Sample_ID'))
delta_abundance[is.na(delta_abundance)] = 0
delta_abundance = delta_abundance %>% mutate(relative_abundance = relative_abundance.y - relative_abundance.x) %>% select(-relative_abundance.x,-relative_abundance.y) %>% mutate(timepoint = 'delta')

abundance_data = bind_rows(abundance_data,delta_abundance)

# get list of all terms to regress over
amr_terms = unique(abundance_data$ARO.Term)

# load and merge metadata

metadata = readRDS('pds08_metadata.rds') #%>% filter(b_bm_weekly<=4.2)
#demarcate responder vs nonresponder
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(rx != 'Treatment',"PLACEBO",if_else(d_bm_weekly>=1,'RESPONDER',if_else(abs(d_bm_weekly)<1,'NOCHANGE-TREATMENT','NOCHANGE-TREATMENT')))))
metadata = metadata %>% mutate(RESP_STATUS = as.factor(if_else(RESP_STATUS != 'PLACEBO',as.character(RESP_STATUS),if_else(d_bm_weekly>=1,'IMPROVED-PLACEBO',if_else(abs(d_bm_weekly)<1,'NOCHANGE-PLACEBO','NOCHANGE-PLACEBO')))))
metadata = metadata %>% mutate(IMPROVED = if_else(RESP_STATUS == 'IMPROVED-PLACEBO' | RESP_STATUS == 'RESPONDER','IMPROVED',if_else(RESP_STATUS == 'NOCHANGE-PLACEBO' | RESP_STATUS == 'NOCHANGE-TREATMENT','NOCHANGE','NOCHANGE')))

metadata = metadata %>% mutate(RESP_V_NONRESP = if_else(RESP_STATUS == 'RESPONDER',1,if_else(RESP_STATUS == 'NOCHANGE-TREATMENT',0,-1))) 
metadata$RESP_V_NONRESP[metadata$RESP_V_NONRESP==-1] = NA

#metadata = metadata %>% filter(b_bm_weekly<5)

# create spoofed metadata table
mdat_amr = list()
for(a in amr_terms){
  b = metadata %>% select(Sample_ID) %>% mutate(ARO.Term = a,timepoint = 'baseline')
  e = metadata %>% select(Sample_ID) %>% mutate(ARO.Term = a,timepoint = 'endpoint')
  d = metadata %>% select(Sample_ID) %>% mutate(ARO.Term = a,timepoint = 'delta')
  samples = bind_rows(b,e,d)
  mdat_amr[[a]] = samples
}
mdat_amr = bind_rows(mdat_amr)

# merge metadata
merged_data = left_join(mdat_amr,metadata,by='Sample_ID')
merged_data = left_join(merged_data,abundance_data,by = c('Sample_ID','ARO.Term','timepoint'))
merged_data$relative_abundance[is.na(merged_data$relative_abundance)] = 0

# for each term, compute 1) endpoint association with rx/resp 2) delta association with rx/resp 
endpoint_rx = list()
endpoint_resp = list()
delta_rx = list()
delta_resp = list()
for(a in amr_terms){
  print(a)
  data_end = merged_data %>% filter(timepoint == 'endpoint',ARO.Term==a)
  if(sum(data_end$relative_abundance)!=0){
    prevalence_rx = (data_end %>% filter(relative_abundance!=0) %>% nrow)/nrow(data_end)
    rx_end = glm(data = data_end, family = 'gaussian',log(relative_abundance + 0.00001) ~ age + b_bm_weekly + rx) %>% tidy %>% filter(term != 'age', term !='b_bm_weekly', term != '(Intercept)') %>% mutate(ARO.Term = a,prevalence = prevalence_rx)
    endpoint_rx[[a]] = rx_end
    prevalence_rnr = (data_end %>% filter(!is.na(RESP_V_NONRESP),relative_abundance!=0) %>% nrow)/nrow(data_end %>% filter(!is.na(RESP_V_NONRESP)))
    if(prevalence_rnr!=0){
      resp_end = glm(data = data_end, family = 'gaussian',log(relative_abundance + 0.00001) ~ age + b_bm_weekly + RESP_V_NONRESP) %>% tidy %>% filter(term != 'age', term !='b_bm_weekly', term != '(Intercept)') %>% mutate(ARO.Term = a,prevalence = prevalence_rnr)
      endpoint_resp[[a]] = resp_end
    }
  }
  data_delta = merged_data %>% filter(timepoint == 'delta',ARO.Term==a)
  if(sum(data_delta$relative_abundance)!=0){
    prevalence_rx = (data_delta %>% filter(relative_abundance!=0) %>% nrow)/nrow(data_delta)
    rx_delta = glm(data = data_delta, family = 'gaussian',log(relative_abundance + abs(min(relative_abundance)) + 0.00001) ~ age + b_bm_weekly + rx) %>% tidy %>% filter(term != 'age', term !='b_bm_weekly', term != '(Intercept)') %>% mutate(ARO.Term = a,prevalence = prevalence_rx)
    delta_rx[[a]] = rx_delta
    prevalence_rnr = (data_delta %>% filter(!is.na(RESP_V_NONRESP),relative_abundance!=0) %>% nrow)/nrow(data_delta %>% filter(!is.na(RESP_V_NONRESP)))
    if(prevalence_rnr!=0){
      resp_delta  = glm(data = data_delta, family = 'gaussian',log(relative_abundance+ abs(min(relative_abundance)) + 0.00001) ~ age + b_bm_weekly + RESP_V_NONRESP) %>% tidy %>% filter(term != 'age', term !='b_bm_weekly', term != '(Intercept)') %>% mutate(ARO.Term = a,prevalence = prevalence_rnr)
      delta_resp[[a]] = resp_delta
    }
  }
}

endpoint_rx_merged = bind_rows(endpoint_rx) %>% mutate(adj = p.adjust(p.value,method = 'BY')) %>% filter(adj<0.1)
endpoint_resp_merged = bind_rows(endpoint_resp) %>% mutate(adj = p.adjust(p.value,method = 'BY')) %>% filter(adj<0.1)
delta_rx_merged = bind_rows(delta_rx) %>% mutate(adj = p.adjust(p.value,method = 'BY')) %>% filter(adj<0.1)
delta_resp_merged = bind_rows(delta_resp) %>% mutate(adj = p.adjust(p.value,method = 'BY')) %>% filter(adj<0.1)

