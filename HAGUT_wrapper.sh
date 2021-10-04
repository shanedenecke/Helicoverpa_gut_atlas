cd /mnt/c/Users/shane/Dropbox/omics_projects/Ha_midgut_atlas
mkdir outputs
mkdir outputs/gene_groups

################################  Protome gather ################################
./scripts/HAGUT_proteome_gather.R

################################  Clean All Omics Data ################################
./scripts/HAGUT_Midgut_atlas_omics_clean.R

################################  Compare Proteome and Transcriptome ################################
./scripts/HAGUT_Proteome_Transcriptome_compare.R

################################  Proteome Venn Diagrams ################################
./scripts/HAGUT_Prot_Venn.R

################################  Transcriptomic groupings ################################
./scripts/HAGUT_Life_stage_DE.R
./scripts/HAGUT_Fuzzy_C_means.R

################################  Protome gather ################################
./scripts/HAGUT_GO_Enrichment.R

################################  Do specific analysis ################################
./scripts/HAGUT_pH_Genes.R
./scripts/HAGUT_P450_Analysis.R
./scripts/HAGUT_SLC_tree.R
