# Immune Cell Subtypes Specific Variants-to-Genes Mapping in IBD

Step 1: Preparing Fine-mapping output
Using fine-mapping script, two fine-mapping datasets of IBD-associated risk loci in EUR were integrated to create a uniform dataset of credible interval SNPs. 

Step 2: Finding SNPs that Impact Binding of TFs
Using atSNP package and JASPAR2024, the EffectSNP script analyzes the impact of credible interval SNPs on the binding affinities of human transcription factors and provides SNP-TF pairs along with their log_lik_ratio. 

Step 3: Finding SNPs that Fall Within Peaks in snATAC-seq
Using snATAC_seq_analysis script, peak-overlapping SNPs could be identified in a cell subtype specific manner. Additionally, This script links risk-mediating peaks to proximal or distal genes through utilizing Cicero single-cell co-accessibility data. 

Step 4: Verifying the Physical Proximity of Co-accessible Pairs of Peaks
Using HiC_enh_evaluation script, the physical proximity of co-accessible peaks from the Cicero single-cell co-accessibility data were screened. 

Step 5: eQTL for Peak-overllaping SNPs
Initially, by using the eQTL_gene_linkage script in eQTL directory, the input for eQTL were prepared. Then using the bash scripts, eQTL was perfoemd through a cell subtype specific manner to find genes associated with SNPs in peaks. 

Step 6: Adding ChromHMM Annotation Labels
ChromHMM_annotation script uses cell subtype specific ChromHMM annotation files to assign the labels of histone marks to risk-mediating peaks. 

Step 7: 
