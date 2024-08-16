#!/bin/bash
#SBATCH --job-name=inspect_seurat
#SBATCH --output=inspect_seurat.out
#SBATCH --error=inspect_seurat.err
#SBATCH --time=010:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load r/4.3.1
module load r-bundle-bioconductor/3.18
Rscript /home/mderakh3/moho/inspect_seurat.R

