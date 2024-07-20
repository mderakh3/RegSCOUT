# Immune Cell Subtypes Specific Variants-to-Genes Mapping in IBD

Table of Contents 
Step 1: Preparing Fine-mapping Output
Step 2: Finding SNPs that Impact Binding of TFs
Step 3: Finding Effect-SNPs that Fall Within Peaks in snATAC-seq
Step 4: Verifying the Physical Proximity of Co-accessible Pairs of Peaks
Step 5: eQTL for Peak-overlapping SNPs
Step 6: Adding ChromHMM Annotation Labels
Step 7: Integration of All Data
Step 8: Prioritizing Maps
Step 9: Cell Subtype Specific Network Analysis
Step 10: Pathway and Disease Enrichment Analysis Using EnrichR
Step 11: Postprocessing
Abbreviations


Step 1: Preparing Fine-mapping Output
Using the fine-mapping script, two fine-mapping datasets of IBD-associated risk loci in EUR individuals were integrated to create a uniform dataset of CI-SNPs.

Step 2: Finding SNPs that Impact Binding of TFs
Using the atSNP package and JASPAR2024, the EffectSNP script analyzes the impact of CI-SNPs on the binding affinities of human TFs and provides EffectSNP-TF pairs along with their log_lik_ratio.

Step 3: Finding Effect-SNPs that Fall Within Peaks in snATAC-seq
Using the snATAC_seq_analysis script, peak-overlapping SNPs can be identified in a cell subtype-specific manner among Effect-SNPs. This script also links risk-mediating peaks to proximal or distal genes using Cicero single-cell co-accessibility data.

Step 4: Verifying the Physical Proximity of Co-accessible Pairs of Peaks
Using the HiC_enh_evaluation script, the physical proximity of co-accessible pairs of peaks from Cicero data, one of which has at least one PO-SNP, were screened.

Step 5: eQTL for Peak-overlapping SNPs
Initially, by using the eQTL_gene_linkage script in the eQTL directory, the input for eQTL containing the list of PO-SNPs was prepared.
Then, using bash scripts, eQTL was performed in a cell subtype-specific manner to find genes associated with SNPs in peaks.

Step 6: Adding ChromHMM Annotation Labels
The ChromHMM_annotation script uses cell subtype-specific ChromHMM annotation files to assign histone mark labels to RMPs.

Step 7: Integration of All Data
Using the final_output script, the data of all variants, their corresponding TFs, RMPs, and genes were integrated in a cell subtype-specific manner.

Step 8: Prioritizing Maps
The prioritizing_maps script uses the final_output and prioritizes those variants-to-genes maps that have genes showing up in at least two methods of gene linkage, and filters out others that do not meet this criterion.

Step 9: Cell Subtype Specific Network Analysis
Using the network_analysis script:
- TFs and genes in each immune cell subtype were isolated to create input for building networks in the STRING database.
- The outputs of STRING were then processed using the same script.

Step 10: Pathway and Disease Enrichment Analysis Using EnrichR
The gene_set_enrichment script performs:
- Pathway analysis using MSigDB.
- Disease analysis using DisGeNet.

Step 11: Postprocessing
In the postprocessing directory:
- part_1, part_2, and part_3 scripts generate tables and figures for the manuscript.

Abbreviations
- SNP: Single-nucleotide polymorphism
- EUR: European
- TF: Transcription factor
- IBD: Inflammatory bowel disease
- CI-SNP: Credible interval SNP
- PO-SNP: Peak-overlapping SNP
- RMP: Risk-mediating peak
- eQTL: Expression quantitative trait loci
- ChromHMM: Chromatin hidden Markov model
