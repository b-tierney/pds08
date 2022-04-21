# Functional response to microbial therapy in the gastrointestinal system of children: a randomized clinical trial

Scripts used in the anaylsis of a paper looking at the impact of a synbiotic on functional constipation in children.

Background
Oral microbial therapy has been studied as an intervention for a range of gastrointestinal disorders. Though research suggests microbial exposure may affect the gastrointestinal system, motility, and host immunity in a pediatric population, data has been inconsistent, with most prior studies being in neither a randomized nor placebo-controlled setting. The aim of this placebo-controlled study was to evaluate the efficacy of a synbiotic on increasing Weekly Bowel Movements (WBMs) in constipated children.

Methods
Sixty-four children (3-17 years of age) were randomized to receive a synbiotic (n=33) comprised of mixed-chain length oligosaccharides and nine microbial strains, or placebo (n=31) for 84 days. Stool microbiota was analyzed on samples collected at baseline and completion. The primary outcome was change from baseline of WBMs in the treatment group compared to placebo.

Results
Treatment increased (p < 0.05) the number of WBMs in children with low WBMs, despite broadly distinctive baseline microbiome signatures. Sequencing revealed that low baseline microbial richness in the treatment group significantly anticipated improvements in constipation (p = 0.00074).

Conclusions
These findings suggest the potential for (i) multi-species-synbiotic interventions to improve digestive health in a pediatric population and (ii) bioinformatics-based methods to predict response to microbial interventions in children.    

### Script Descriptions

#### build_metadata_and_process_abundances.R
Parse abundance data from MetaPhlAn3/HUMAnN3 and merge clinical + sequencing metadata.

#### metaphlan_multilevel_analysis.R
Compute and plot richness at different taxonomic levels across different samples between baseline and endpoint.

#### bracken_multi_level_analysis.R
Compute and plot richness at different taxonomic levels across different samples between baseline and endpoint using an alternative taxonomic framework (Kraken2/Bracken).

#### metaphlan_shannon_multilevel.R
Compute shannon diversity at different taxonomic levels based on MetaPhlAn3 output. Not used in manuscript.

#### metaphlan_simpson_multilevel.R
Compute simpson diversity at different taxonomic levels based on MetaPhlAn3 output. Not used in manuscript.

#### strain_detection.R
Merge output of Anvio-7 to identify abundances of PDS-08 genomes in each sample.

#### bifidogenic.R
Merge output of Anvio-7 to identify abundances of Bifidobacterial genomes in each sample.

#### taxonomic_regressions_clinical_data.R
Compute associations with different clinical variables and HUMAnN3 pathway abundances.

#### pathway_regressions_clinical_data.R
Compute associations with different clinical variables and HUMAnN3 pathway abundances.
