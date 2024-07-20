# Immune Cell Subtypes Specific Variants-to-Genes Mapping in IBD

Step 1: Preparing Fine-mapping output
Using fine-mapping script, two fine-mapping datasets of IBD-associated risk loci in EUR individuals were integrated to create a uniform dataset of CI-SNPs. 

Step 2: Finding SNPs that Impact Binding of TFs
Using atSNP package and JASPAR2024, the EffectSNP script analyzes the impact of CI-SNPs on the binding affinities of human TFs and provides EffectSNP-TF pairs along with their log_lik_ratio. 

Step 3: Finding Effect-SNPs that Fall Within Peaks in snATAC-seq
Using snATAC_seq_analysis script, peak-overlapping SNPs could be identified in a cell subtype specific manner among Effect-SNPs. Additionally, this script links risk-mediating peaks to proximal or distal genes through utilizing Cicero single-cell co-accessibility data. 

Step 4: Verifying the Physical Proximity of Co-accessible Pairs of Peaks
Using HiC_enh_evaluation script, the physical proximity of co-accessible pairs of peaks from Cicero data, one of which has at least on PO-SNPs were screened. 

Step 5: eQTL for Peak-overllaping SNPs
Initially, by using the eQTL_gene_linkage script in eQTL directory, the input for eQTL which contains the list PO-SNPs were prepared. Then using the bash scripts, eQTL was perfoemd through a cell subtype specific manner to find genes associated with SNPs in peaks. 

Step 6: Adding ChromHMM Annotation Labels
ChromHMM_annotation script uses cell subtype specific ChromHMM annotation files to assign the labels of histone marks to RMPs. 

Step 7: Integration of All Data
Using the final_output script, the data of all variants, their corresponding TFs, RMPs, and genes were integrated in a cell subtype specific manner. 

Step 8: Prioritizing Maps
The prioritizing_maps script uses the final_ouput output and prioritizes those variants-to-genes maps that have genes that showed up in two methods of gene linkage at least and it gets rid of other genes and those maps that do not pass this criteria. 

Step 9: Cell Subtype Specific Network Analysis
Using the network_analysis script, TFs and genes in each immune cell subtype were isolated to create input for building networks in STRING database. Then, using the same script, the outputs of STRING were processed. 

Step 10: Pathway and Disease Enrichement Analysis Using EnrichR
The gene_set_enrichment script uses the all TFs and genes in each immune cell subtype separately and performs pathway analysis using MSigDB and disease analysis using DisGeNet. 

Step 11: Postprocessing
In the postprocessing directory, part_1, part_2, and part_3 scripts generate tables and figures for the manuscript. 

Abbreviations:
SNP; Single-nucleotide polymorphism
EUR; European
TF: Transcription factor
IBD; Inflammatory bowel disease
CI-SNP; Credible interval SNP
PO-SNP; Peak-overlapping SNP
RMP; Risk-mediating peak
eQTL; Expression quantitative trait loci
ChromHMM; Chromatin hidden Markov model

