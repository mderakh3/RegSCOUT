#!/bin/bash
#SBATCH --job-name=filter_seurat
#SBATCH --output=filter_seurat.out
#SBATCH --error=filter_seurat.err
#SBATCH --time=010:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load r/4.3.1
module load r-bundle-bioconductor/3.18
Rscript /home/mderakh3/moho/filter_seurat.R
