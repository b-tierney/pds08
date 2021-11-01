#!/usr/bin/Rscript --vanilla

# generate plots for regression output

library(tidyverse)
library(ggplot2)
library(umap)
library(cowplot)
library(vegan)
library(ggrepel)

theme_set(theme_cowplot())

args <- commandArgs(trailingOnly = TRUE)

system('mkdir metaphlan_associations/volcano_plots_baseline')
system('mkdir metaphlan_associations/volcano_plots_delta')
system('mkdir metaphlan_associations/volcano_plots_endpoint')

system('mkdir metaphlan_associations/volcano_plots_baseline_interaction')
system('mkdir metaphlan_associations/volcano_plots_delta_interaction')
system('mkdir metaphlan_associations/volcano_plots_endpoint_interaction')

### alpha diversity

baseline = list.files('metaphlan_associations')
baseline = baseline[grep('baseline',baseline)]
baseline = baseline[grep('diversity',baseline)]
baseline = baseline[grep('metaphlan',baseline)]
baseline = baseline[grep('outcome',baseline)]

data = list()
for(f in baseline){
  d = readRDS(paste('metaphlan_associations/',f,sep='')) %>% mutate(var = gsub('_metaphlan_baseline_diversity.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

base = bind_rows(data)

delta = list.files('metaphlan_associations')
delta = delta[grep('delta',delta)]
delta = delta[grep('diversity',delta)]
delta = delta[grep('metaphlan',delta)]
delta = delta[grep('outcome',delta)]

data = list()
for(f in delta){
  d = readRDS(paste('metaphlan_associations/',f,sep='')) %>% mutate(var = gsub('_metaphlan_delta_diversity.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

delta = bind_rows(data)

endpoint = list.files('metaphlan_associations')
endpoint = endpoint[grep('endpoint',endpoint)]
endpoint = endpoint[grep('diversity',endpoint)]
endpoint = endpoint[grep('metaphlan',endpoint)]
endpoint = endpoint[grep('outcome',endpoint)]


data = list()
for(f in endpoint){
  d = readRDS(paste('metaphlan_associations/',f,sep='')) %>% mutate(var = gsub('_metaphlan_endpoint_diversity.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

end = bind_rows(data)

alpha = bind_rows(base%>% mutate(tp = 'baseline'),delta%>% mutate(tp = 'delta'),end %>% mutate(tp = 'endpoint'))
write.csv(alpha,'metaphlan_associations/diversity_clinical.csv')

baseline_alpha = readRDS('./metaphlan_associations/regression_output_microbe_treatment_metaphlan_baseline_diversity.rds') %>% mutate()
endpoint_alpha = readRDS('./metaphlan_associations/regression_output_microbe_treatment_metaphlan_endpoint_diversity.rds') %>% mutate()
delta_alpha = readRDS('./metaphlan_associations/regression_output_microbe_treatment_metaphlan_delta_diversity.rds') %>% mutate()

treat = bind_rows(baseline_alpha %>% mutate(tp = 'baseline'),delta_alpha %>% mutate(tp = 'delta'),endpoint_alpha %>% mutate(tp = 'endpoint'))
write.csv(treat,'metaphlan_associations/diversity_treatment.csv')

files = list.files('metaphlan_associations')
files = files[grep('baseline',files)]
files = files[grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste('metaphlan_associations/',f,sep='')) %>% mutate(var = gsub('_metaphlan_endpoint_diversity.rds','',gsub('regression_output_interaction_microbe_','',f)))
  data[[f]] =d
}

base = bind_rows(data)

files = list.files('metaphlan_associations')
files = files[grep('endpoint',files)]
files = files[grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste('metaphlan_associations/',f,sep='')) %>% mutate(var = gsub('_metaphlan_endpoint_diversity.rds','',gsub('regression_output_interaction_microbe_','',f)))
  data[[f]] =d
}

end = bind_rows(data)

files = list.files('metaphlan_associations')
files = files[grep('delta',files)]
files = files[grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste('metaphlan_associations/',f,sep='')) %>% mutate(var = gsub('_metaphlan_endpoint_diversity.rds','',gsub('regression_output_interaction_microbe_','',f)))
  data[[f]] =d
}

delta = bind_rows(data)

alpha_int = bind_rows(base%>% mutate(tp = 'baseline'),delta%>% mutate(tp = 'delta'),end %>% mutate(tp = 'endpoint'))
write.csv(alpha_int,'metaphlan_associations/diversity_clinical_interaction.csv')

get_volcanos <- function(regression_output,output_folder){
  features = unique(regression_output$var)
  for(f in features){
    data_sub = regression_output %>% filter(var==f)
    plotout = ggplot(data=data_sub,aes(x=estimate,y=-log10(bh),label=term))+geom_point()+geom_hline(yintercept = -log10(0.1),linetype='dashed') + ggtitle(f) + ylab('-log10(qval)') +ggtitle('') + xlim(-5,5) + ylim(0,11) + geom_label_repel(data = data_sub %>% arrange(bh) %>% filter(bh<0.1) %>% head(10),aes(label = term),
                  box.padding   = 0.1, 
                  point.padding = 0.1,
                  max.overlaps=100,
                  size=3,
                  segment.color = 'grey50')
    ggsave(paste(output_folder,f,'_volcano.pdf',sep=''),width=10,height=10)
  }
}

setwd('metaphlan_associations')

files = list.files()
files = files[grep('baseline',files)]
files = files[-grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('outcome',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_metaphlan_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

baseline_metaphlan_merged = bind_rows(data)
baseline_metaphlan_merged = baseline_metaphlan_merged %>% filter(term!='k__Bacteria')
baseline_metaphlan_merged = baseline_metaphlan_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(baseline_metaphlan_merged,'volcano_plots_baseline/')

files = list.files()
files = files[grep('endpoint',files)]
files = files[-grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('outcome',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_metaphlan_endpoint.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

endpoint_metaphlan_merged = bind_rows(data)
endpoint_metaphlan_merged = endpoint_metaphlan_merged %>% filter(term!='k__Bacteria')
endpoint_metaphlan_merged = endpoint_metaphlan_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(endpoint_metaphlan_merged,'volcano_plots_endpoint/')

files = list.files()
files = files[grep('delta',files)]
files = files[-grep('diversity',files)]
files = files[-grep('plots',files)]
files = files[grep('metaphlan',files)]
files = files[grep('outcome',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_metaphlan_delta.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

delta_metaphlan_merged = bind_rows(data)
delta_metaphlan_merged = delta_metaphlan_merged %>% filter(term!='k__Bacteria')
delta_metaphlan_merged = delta_metaphlan_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(delta_metaphlan_merged,'volcano_plots_delta/')

output = bind_rows(baseline_metaphlan_merged %>% mutate(tp = 'baseline'),delta_metaphlan_merged%>% mutate(tp = 'delta'),endpoint_metaphlan_merged%>% mutate(tp = 'endpoint'))
write.csv(output,'regression_output.csv')

### INTERACTION MODEL

files = list.files()
files = files[grep('baseline',files)]
files = files[-grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_metaphlan_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

baseline_metaphlan_merged = bind_rows(data)
baseline_metaphlan_merged = baseline_metaphlan_merged %>% filter(term!='k__Bacteria')
baseline_metaphlan_merged = baseline_metaphlan_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(baseline_metaphlan_merged,'volcano_plots_baseline_interaction/')

files = list.files()
files = files[grep('endpoint',files)]
files = files[-grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_metaphlan_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

endpoint_metaphlan_merged = bind_rows(data)
endpoint_metaphlan_merged = endpoint_metaphlan_merged %>% filter(term!='k__Bacteria')
endpoint_metaphlan_merged = endpoint_metaphlan_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(endpoint_metaphlan_merged,'volcano_plots_endpoint_interaction/')

files = list.files()
files = files[grep('delta',files)]
files = files[-grep('diversity',files)]
files = files[grep('metaphlan',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_metaphlan_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

delta_metaphlan_merged = bind_rows(data)
delta_metaphlan_merged = delta_metaphlan_merged %>% filter(term!='k__Bacteria')
delta_metaphlan_merged = delta_metaphlan_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(delta_metaphlan_merged,'volcano_plots_delta_interaction/')

### treatment

baseline = readRDS('regression_output_microbe_treatment_metaphlan_baseline.rds') 
endpoint = readRDS('regression_output_microbe_treatment_metaphlan_endpoint.rds')
delta = readRDS('regression_output_microbe_treatment_metaphlan_delta.rds') 

treatment_regression = bind_rows(baseline%>% mutate(tp = 'baseline'),delta%>% mutate(tp = 'delta'),endpoint %>% mutate(tp = 'endpoint'))
write.csv(treatment_regression,'./treatment_regressions.csv')

###### PATHWAYS
setwd('../')

system('mkdir pathways_associations/volcano_plots_baseline')
system('mkdir pathways_associations/volcano_plots_delta')
system('mkdir pathways_associations/volcano_plots_endpoint')

system('mkdir pathways_associations/volcano_plots_baseline_interaction')
system('mkdir pathways_associations/volcano_plots_delta_interaction')
system('mkdir pathways_associations/volcano_plots_endpoint_interaction')

### alpha diversity

baseline = list.files('pathways_associations')
baseline = baseline[grep('baseline',baseline)]
baseline = baseline[grep('diversity',baseline)]
baseline = baseline[grep('pathways',baseline)]
baseline = baseline[grep('outcome',baseline)]

data = list()
for(f in baseline){
  d = readRDS(paste('pathways_associations/',f,sep='')) %>% mutate(var = gsub('_pathways_baseline_diversity.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

base = bind_rows(data)

delta = list.files('pathways_associations')
delta = delta[grep('delta',delta)]
delta = delta[grep('diversity',delta)]
delta = delta[grep('pathways',delta)]
delta = delta[grep('outcome',delta)]

data = list()
for(f in delta){
  d = readRDS(paste('pathways_associations/',f,sep='')) %>% mutate(var = gsub('_pathways_delta_diversity.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

delta = bind_rows(data)

endpoint = list.files('pathways_associations')
endpoint = endpoint[grep('endpoint',endpoint)]
endpoint = endpoint[grep('diversity',endpoint)]
endpoint = endpoint[grep('pathways',endpoint)]
endpoint = endpoint[grep('outcome',endpoint)]


data = list()
for(f in endpoint){
  d = readRDS(paste('pathways_associations/',f,sep='')) %>% mutate(var = gsub('_pathways_endpoint_diversity.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

end = bind_rows(data)

alpha = bind_rows(base%>% mutate(tp = 'baseline'),delta%>% mutate(tp = 'delta'),end %>% mutate(tp = 'endpoint'))
write.csv(alpha,'pathways_associations/diversity_clinical.csv')

baseline_alpha = readRDS('./pathways_associations/regression_output_microbe_treatment_pathways_baseline_diversity.rds') %>% mutate()
endpoint_alpha = readRDS('./pathways_associations/regression_output_microbe_treatment_pathways_endpoint_diversity.rds') %>% mutate()
delta_alpha = readRDS('./pathways_associations/regression_output_microbe_treatment_pathways_delta_diversity.rds') %>% mutate()

treat = bind_rows(baseline_alpha %>% mutate(tp = 'baseline'),delta_alpha %>% mutate(tp = 'delta'),endpoint_alpha %>% mutate(tp = 'endpoint'))
write.csv(treat,'pathways_associations/diversity_treatment.csv')

files = list.files('pathways_associations')
files = files[grep('baseline',files)]
files = files[grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste('pathways_associations/',f,sep='')) %>% mutate(var = gsub('_pathways_endpoint_diversity.rds','',gsub('regression_output_interaction_microbe_','',f)))
  data[[f]] =d
}

base = bind_rows(data)

files = list.files('pathways_associations')
files = files[grep('endpoint',files)]
files = files[grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste('pathways_associations/',f,sep='')) %>% mutate(var = gsub('_pathways_endpoint_diversity.rds','',gsub('regression_output_interaction_microbe_','',f)))
  data[[f]] =d
}

end = bind_rows(data)

files = list.files('pathways_associations')
files = files[grep('delta',files)]
files = files[grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste('pathways_associations/',f,sep='')) %>% mutate(var = gsub('_pathways_endpoint_diversity.rds','',gsub('regression_output_interaction_microbe_','',f)))
  data[[f]] =d
}

delta = bind_rows(data)

alpha_int = bind_rows(base%>% mutate(tp = 'baseline'),delta%>% mutate(tp = 'delta'),end %>% mutate(tp = 'endpoint'))
write.csv(alpha_int,'pathways_associations/diversity_clinical_interaction.csv')

get_volcanos <- function(regression_output,output_folder){
  features = unique(regression_output$var)
  for(f in features){
    data_sub = regression_output %>% filter(var==f)
    plotout = ggplot(data=data_sub,aes(x=estimate,y=-log10(bh),label=term))+geom_point()+geom_hline(yintercept = -log10(0.05),linetype='dashed') + ggtitle(f) + ylab('-log10(qval)') +ggtitle('') + xlim(-5,5) + ylim(0,11) + geom_label_repel(data = data_sub %>% arrange(bh) %>% filter(bh<0.1) %>% head(10),aes(label = term),
                  box.padding   = 0.1, 
                  point.padding = 0.1,
                  max.overlaps=100,
                  size=3,
                  segment.color = 'grey50')
    ggsave(paste(output_folder,f,'_volcano.pdf',sep=''),width=10,height=10)
  }
}

setwd('pathways_associations')

files = list.files()
files = files[grep('baseline',files)]
files = files[-grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('outcome',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_pathways_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

baseline_pathways_merged = bind_rows(data)
baseline_pathways_merged = baseline_pathways_merged %>% filter(term!='k__Bacteria')
baseline_pathways_merged = baseline_pathways_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(baseline_pathways_merged,'volcano_plots_baseline/')

files = list.files()
files = files[grep('endpoint',files)]
files = files[-grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('outcome',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_pathways_endpoint.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

endpoint_pathways_merged = bind_rows(data)
endpoint_pathways_merged = endpoint_pathways_merged %>% filter(term!='k__Bacteria')
endpoint_pathways_merged = endpoint_pathways_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(endpoint_pathways_merged,'volcano_plots_endpoint/')

files = list.files()
files = files[grep('delta',files)]
files = files[-grep('diversity',files)]
files = files[-grep('plots',files)]
files = files[grep('pathways',files)]
files = files[grep('outcome',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_pathways_delta.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

delta_pathways_merged = bind_rows(data)
delta_pathways_merged = delta_pathways_merged %>% filter(term!='k__Bacteria')
delta_pathways_merged = delta_pathways_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(delta_pathways_merged,'volcano_plots_delta/')

output = bind_rows(baseline_pathways_merged %>% mutate(tp = 'baseline'),delta_pathways_merged%>% mutate(tp = 'delta'),endpoint_pathways_merged%>% mutate(tp = 'endpoint'))
write.csv(output,'regression_output.csv')

### INTERACTION MODEL

files = list.files()
files = files[grep('baseline',files)]
files = files[-grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_pathways_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

baseline_pathways_merged = bind_rows(data)
baseline_pathways_merged = baseline_pathways_merged %>% filter(term!='k__Bacteria')
baseline_pathways_merged = baseline_pathways_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(baseline_pathways_merged,'volcano_plots_baseline_interaction/')

files = list.files()
files = files[grep('endpoint',files)]
files = files[-grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_pathways_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

endpoint_pathways_merged = bind_rows(data)
endpoint_pathways_merged = endpoint_pathways_merged %>% filter(term!='k__Bacteria')
endpoint_pathways_merged = endpoint_pathways_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(endpoint_pathways_merged,'volcano_plots_endpoint_interaction/')

files = list.files()
files = files[grep('delta',files)]
files = files[-grep('diversity',files)]
files = files[grep('pathways',files)]
files = files[grep('interaction',files)]

data = list()
for(f in files){
  d = readRDS(paste(f,sep='')) %>% mutate(var = gsub('_pathways_baseline.rds','',gsub('regression_output_outcome_microbe_','',f)))
  data[[f]] =d
}

delta_pathways_merged = bind_rows(data)
delta_pathways_merged = delta_pathways_merged %>% filter(term!='k__Bacteria')
delta_pathways_merged = delta_pathways_merged %>% mutate(term = gsub("^.*\\|", "", term))

get_volcanos(delta_pathways_merged,'volcano_plots_delta_interaction/')

### treatment

baseline = readRDS('regression_output_microbe_treatment_pathways_baseline.rds') 
endpoint = readRDS('regression_output_microbe_treatment_pathways_endpoint.rds')
delta = readRDS('regression_output_microbe_treatment_pathways_delta.rds') 

treatment_regression = bind_rows(baseline%>% mutate(tp = 'baseline'),delta%>% mutate(tp = 'delta'),endpoint %>% mutate(tp = 'endpoint'))
write.csv(treatment_regression,'./treatment_regressions.csv')

