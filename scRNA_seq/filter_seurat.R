#!/usr/bin/env Rscript

library(Seurat)

# Load the Seurat object
seurat_object <- readRDS('/home/mderakh3/moho/scIBD.gex_matrix.rds')

# Define the samples to keep based on your screenshot
samples_to_keep <- c("GSM4761136", "GSM4761137", "GSM4761138", "GSM4761139", 
                     "GSM4761140", "GSM4761141", "GSM4761142", "GSM4761143", 
                     "GSM4761144", "GSM3576411", "GSM3576412", "GSM3576413", 
                     "GSM3576414", "GSM3576415", "GSM3576416", "GSM3576417", 
                     "GSM3576418", "GSM3576419", "GSM3576420", "GSM3576421", 
                     "GSM3576422", "GSM3576423", "GSM3576424", "GSM3576425")

# Filter the Seurat object to keep only the cells belonging to these samples
cells_to_keep <- WhichCells(seurat_object, expression = sample %in% samples_to_keep)

# Subset the Seurat object
filtered_seurat_object <- subset(seurat_object, cells = cells_to_keep)

# Save the filtered Seurat object
saveRDS(filtered_seurat_object, file = '/home/mderakh3/moho/scIBD_filtered.rds')
