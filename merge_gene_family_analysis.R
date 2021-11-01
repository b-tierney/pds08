# merge gene family analysis

library(tidyverse)
library(stringi)

files=list.files()

basenames  = map(files, function(x) gsub('regression_output_outcome_microbe_','',x) %>% gsub("_genefamilies_baseline.rds","",.) %>% strsplit(., "") %>% unlist %>% rev(.) %>% paste(.,collapse="") %>% gsub("[^_]*_(.*)", "\\1", .) %>% strsplit(., "") %>% unlist %>% rev(.) %>% paste(.,collapse="")) %>% unlist %>% unname 

mapping = tibble(files,basenames)

for(b in basenames){
	print(b)
	subfiles = mapping %>% filter(basenames==b) %>% select(files)
	data = list()
	for(f in files){
		data[[f]] = readRDS(f)
	}
	data = bind_rows(data) %>% mutate(bh = p.adjust(p.value,method ='BH')) 
	output[[b]] = data
}

output = bind_rows(output)







