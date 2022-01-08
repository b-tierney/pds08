# Functional response to microbial therapy in the gastrointestinal system of children

Scripts used in the anaylsis of a paper looking at the impact of a synbiotic on functional constipation in children.

Recent advances in bioinformatics have enabled a deeper understanding of microbial ecology and the gut microbiome’s role in human health and wellbeing. Gastrointestinal microbes exert functional influence on the host through a range of metabolic and immunological mechanisms and, in exchange, the host shapes resident microbial communities through diet, nutrition, lifestyle, and medication. Oral microbial therapy has been extensively studied from infancy through adulthood with inconsistent and variable results, though emerging research suggests microbial exposure may intimately affect the gastrointestinal system, neuromotility, and host immunity. Here, we report the effects of a defined microbial consortia suspended in a preferential substrate of mixed chain oligosaccharides on the gastrointestinal microbiome and gastric motility in children compared to placebo. Microbial ingestion significantly improved the number of spontaneous bowel movements (BMs) in children irrespective of highly personalized microbiome signatures at baseline. Metagenomic shotgun sequencing revealed that baseline microbial richness was able to distinguish between responders and non-responders. Microbial richness alone significantly anticipated improvements in constipation in responders versus non responders (p<0.001). Overall, this research supports the potential for microbial interventions to improve digestive health in a pediatric population and bioinformatics-based methods to predict personalized responses to microbial interventions in children.

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
