#!/usr/bin/Rscript --vanilla

# compute alpha delta


args <- commandArgs(trailingOnly = TRUE)

base = args[[1]]
end = args[[2]]

base = 'metaphlan_baseline_diversity.rds'
end = 'metaphlan_endpoint_diversity.rds'