#!/bin/bash

# Define the source directory
SOURCE_DIR="/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage/eQTL_input"

# Define the destination
DESTINATION="mderakh3@gra-login1.computecanada.ca:~/moho/"

# List of files to transfer
FILES=(
    "cMono_eqtl_snps.txt"
    "cd16_Mono_eqtl_snps.txt"
    "act_cd4_t_eqtl_snps.txt"
    "cyto_cd8_t_eqtl_snps.txt"
    "mem_b_eqtl_snps.txt"
    "mem_cd8_t_eqtl_snps.txt"
    "naive_b_eqtl_snps.txt"
    "naive_cd4_t_eqtl_snps.txt"
    "naive_cd8_t_eqtl_snps.txt"
    "tReg_eqtl_snps.txt"
    "nk_eqtl_snps.txt"
)

# Loop through the files and transfer each one
for FILE in "${FILES[@]}"; do
    scp "$SOURCE_DIR/$FILE" "$DESTINATION"
done

echo "All files have been transferred."
