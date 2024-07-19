#!/bin/bash

# Define the path to the tabix binary
tabix_path="/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/tabix/0.2.6/bin/tabix"

# Check if the tabix binary exists
if [ ! -x "$tabix_path" ]; then
    echo "Tabix binary not found or not executable at $tabix_path"
    exit 1
fi

# Function to process each SNP file
process_snps() {
    local input_file=$1
    local url=$2
    local output_file=$3

    if [ ! -f "$input_file" ]; then
        echo "Input file $input_file not found"
        return 1
    fi

    xargs -a "$input_file" -I {} $tabix_path -fh "$url" {} > "$output_file"
}

# Process each file
process_snps /home/mderakh3/moho/nk_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000014/QTD000115/QTD000115.all.tsv.gz \
    /home/mderakh3/moho/nk_eqtl.txt

process_snps /home/mderakh3/moho/tReg_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000469/QTD000469.all.tsv.gz \
    /home/mderakh3/moho/tReg_eqtl_mem.txt

process_snps /home/mderakh3/moho/tReg_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000003/QTD000036/QTD000036.all.tsv.gz \
    /home/mderakh3/moho/tReg_eqtl_nai.txt

process_snps /home/mderakh3/moho/naive_cd8_t_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000489/QTD000489.all.tsv.gz \
    /home/mderakh3/moho/naive_cd8_t_eqtl.txt

process_snps /home/mderakh3/moho/naive_cd4_t_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000479/QTD000479.all.tsv.gz \
    /home/mderakh3/moho/naive_cd4_t_eqtl.txt

process_snps /home/mderakh3/moho/nai_b_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000474/QTD000474.all.tsv.gz \
    /home/mderakh3/moho/naive_b_eqtl.txt

process_snps /home/mderakh3/moho/mem_cd8_t_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000007/QTD000066/QTD000066.all.tsv.gz \
    /home/mderakh3/moho/mem_cd8_t_eqtl.txt

process_snps /home/mderakh3/moho/mem_b_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000009/QTD000080/QTD000080.all.tsv.gz \
    /home/mderakh3/moho/mem_b_eqtl.txt

process_snps /home/mderakh3/moho/cyto_cd8_t_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000494/QTD000494.all.tsv.gz \
    /home/mderakh3/moho/cyto_cd8_t_eqtl.txt

process_snps /home/mderakh3/moho/act_cd4_t_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000484/QTD000484.all.tsv.gz \
    /home/mderakh3/moho/act_cd4_t_eqtl.txt

process_snps /home/mderakh3/moho/cd16_Mono_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000499/QTD000499.all.tsv.gz \
    /home/mderakh3/moho/cd16_Mono_eqtl.txt

process_snps /home/mderakh3/moho/cMono_eqtl_snps.txt \
    ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000007/QTD000069/QTD000069.all.tsv.gz \
    /home/mderakh3/moho/cMono_eqtl.txt

echo "All files have been transferred."