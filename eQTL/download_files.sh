#!/bin/bash

# Define the destination directory
DEST_DIR="/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage/eQTL_output"

# Define the source
SOURCE="mderakh3@gra-login1.computecanada.ca:~/moho/"

# List of files to download
FILES=(
    "cMono_eqtl.txt"
    "cd16_Mono_eqtl.txt"
    "act_cd4_t_eqtl.txt"
    "cyto_cd8_t_eqtl.txt"
    "mem_b_eqtl.txt"
    "naive_b_eqtl.txt"
    "naive_cd4_t_eqtl.txt"
    "naive_cd8_t_eqtl.txt"
    "tReg_eqtl_nai.txt"
    "tReg_eqtl_mem.txt"
    "nk_eqtl.txt"
)

# Loop through the files and download each one
for FILE in "${FILES[@]}"; do
    scp "$SOURCE$FILE" "$DEST_DIR/"
done

echo "All files have been downloaded."