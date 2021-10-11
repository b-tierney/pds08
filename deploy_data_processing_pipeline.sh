#!/bin/bash

# run data prep pipeline

Rscript build_metadata_and_process_abundances.R

for file in metaphlan*rds; do Rscript ~/Github/pds08/compute_alpha_and_beta_diversities.R $file; echo $file; done
Rscript compute_alpha_delta.R metaphlan_baseline_diversity.rds metaphlan_endpoint_diversity.rds

for file in pathways*rds; do Rscript ~/Github/pds08/compute_alpha_and_beta_diversities.R $file; echo $file; done
Rscript compute_alpha_delta.R pathways_baseline_diversity.rds pathways_endpoint_diversity.rds

for file in genefamilies*rds; do Rscript ~/Github/pds08/compute_alpha_and_beta_diversities.R $file; echo $file; done
Rscript compute_alpha_delta.R genefamilies_baseline_diversity.rds genefamilies_endpoint_diversity.rds

