#!/bin/bash

# deploy regressions

while read p; do Rscript ~/GitHub/pds08/compute_microbe-treatment_associations.R pds08_metadata.rds $p; done<metaphlan_data

while read pp; do while read p; do Rscript ~/GitHub/pds08/compute_clinical-microbe_associations.R pds08_metadata.rds $p $pp; done<metaphlan_data; done<clinical_variables


