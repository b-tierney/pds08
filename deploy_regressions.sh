#!/bin/bash

# deploy regressions

# metaphlan
echo 'COMPUTING TAXONOMIC ASSOCIATIONS'

echo 'Treatment regressions'
mkdir metaphlan_associations
while read p; do Rscript compute_microbe-treatment_associations.R pds08_metadata.rds $p metaphlan; done<metaphlan_data
echo 'Outcome regressions'
while read pp; do while read p; do echo $pp; Rscript compute_clinical-microbe_associations.R pds08_metadata.rds $p $pp metaphlan; done<metaphlan_data; done<clinical_variables
echo 'Interaction regressions'
while read pp; do while read p; do Rscript compute_interaction_model_associations.R pds08_metadata.rds $p $pp metaphlan; done<metaphlan_data; done<clinical_variables

# pathways

echo 'COMPUTING PATHWAY ASSOCIATIONS'
mkdir pathways_associations
echo 'Treatment regressions'
while read p; do Rscript compute_microbe-treatment_associations.R pds08_metadata.rds $p pathways; done<pathways_data
echo 'Outcome regressions'
while read pp; do while read p; do echo $pp; Rscript compute_clinical-microbe_associations.R pds08_metadata.rds $p $pp pathways; done<pathways_data; done<clinical_variables
echo 'Interaction regressions'
while read pp; do while read p; do Rscript compute_interaction_model_associations.R pds08_metadata.rds $p $pp pathways; done<pathways_data; done<clinical_variables

### NEED TO CHUNK GENE FAMILIES AND THEN POSTPROCESS WITH APPROPRIATE SCRIPT˜

# gene families
echo 'COMPUTING GENE FAMILY ASSOCIATIONS'
mkdir genefamilies_associations 
echo 'Treatment regressions'
while read p; do Rscript compute_microbe-treatment_associations.R pds08_metadata.rds $p genefamilies; done<genefamilies_data
echo 'Outcome regressions'
while read pp; do while read p; do echo $pp; Rscript compute_clinical-microbe_associations.R pds08_metadata.rds $p $pp genefamilies; done<genefamilies_data; done<clinical_variables
echo 'Interaction regressions'
while read pp; do while read p; do Rscript compute_interaction_model_associations.R pds08_metadata.rds $p $pp genefamilies; done<genefamilies_data; done<clinical_variables


