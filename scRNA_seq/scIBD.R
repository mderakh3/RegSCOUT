setwd("/Users/m.hossein_drn/Documents/aim_1/scRNA_seq/")

library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(openxlsx)

# reading the expression data
scIBD_ids <- readxl::read_excel("scIBD_ids.xlsx", sheet = 2) #PBMC samples selected
scIBD_seurat <- readRDS("scIBD_filtered.rds")

# preparing genes from GWAS
cell_subtypes_to_types <- c(
  act_cd4_t = "CD4T",
  naive_cd4_t = "CD4T",
  tReg = "CD4T",
  cyto_cd8_t = "CD8T",
  naive_cd8_t = "CD8T",
  mem_cd8_t = "CD8T",
  mem_b = "B_Plasma",
  naive_b = "B_Plasma",
  iMono = "Myeloid",
  ncMono = "Myeloid",
  cMono = "Myeloid",
  adaptive_NK = "ILC",
  cyto_nk = "ILC",
  mkc = "Myeloid",
  pDC = "Myeloid",
  cDC = "Myeloid",
  plasma = "B_Plasma"
)

scored_genes <- read.delim("scored_table.txt") # all genes

scored_genes <- scored_genes %>% 
  dplyr::select(trait, gene, cell_type) %>%
  distinct()

scored_genes$cluster <- sapply(scored_genes$cell_type, function(x) {
  cell_subtypes_to_types[x]
})

prio_genes <- read.delim("prioritized_table.txt") # prioritized genes
prio_genes <- prio_genes %>% separate_rows(gene_score, sep = ", ")
prio_genes$gene_score <- gsub("\\(.*\\)", "", prio_genes$gene_score)
prio_genes$gene_score <- gsub(" ", "", prio_genes$gene_score)

prio_genes <- prio_genes %>%
  dplyr::select(trait, gene_score, cell_type) %>%
  distinct() %>%
  rename(gene = gene_score)

prio_genes$cluster <- sapply(prio_genes$cell_type, function(x) {
  cell_subtypes_to_types[x]
})

# subset the data
for (i in unique(scIBD_seurat$major_cluster)) {
  assign(paste0(i, "_cells"), subset(scIBD_seurat, major_cluster == i))
}

Idents(Myeloid_cells) <- Myeloid_cells$disease
Idents(CD4T_cells) <- CD4T_cells$disease
Idents(CD8T_cells) <- CD8T_cells$disease
Idents(ILC_cells) <- ILC_cells$disease
Idents(B_Plasma_cells) <- B_Plasma_cells$disease

# finding DEGs
UC_Myeloid_degs <- FindMarkers(Myeloid_cells, ident.1 = "UC_PBMC", ident.2 = "Healthy_PBMC")
CD_Myeloid_degs <- FindMarkers(Myeloid_cells, ident.1 = "CD_PBMC", ident.2 = "Healthy_PBMC")

UC_CD4T_degs <- FindMarkers(CD4T_cells, ident.1 = "UC_PBMC", ident.2 = "Healthy_PBMC")
CD_CD4T_degs <- FindMarkers(CD4T_cells, ident.1 = "CD_PBMC", ident.2 = "Healthy_PBMC")

UC_CD8T_degs <- FindMarkers(CD8T_cells, ident.1 = "UC_PBMC", ident.2 = "Healthy_PBMC")
CD_CD8T_degs <- FindMarkers(CD8T_cells, ident.1 = "CD_PBMC", ident.2 = "Healthy_PBMC")

UC_ILC_degs <- FindMarkers(ILC_cells, ident.1 = "UC_PBMC", ident.2 = "Healthy_PBMC")
CD_ILC_degs <- FindMarkers(ILC_cells, ident.1 = "CD_PBMC", ident.2 = "Healthy_PBMC")

UC_B_Plasma_degs <- FindMarkers(B_Plasma_cells, ident.1 = "UC_PBMC", ident.2 = "Healthy_PBMC")
CD_B_Plasma_degs <- FindMarkers(B_Plasma_cells, ident.1 = "CD_PBMC", ident.2 = "Healthy_PBMC")

# checking for all genes in the DEGs
dataframe_map <- list(
  UC_Myeloid_degs = UC_Myeloid_degs,
  CD_Myeloid_degs = CD_Myeloid_degs,
  UC_CD4T_degs = UC_CD4T_degs,
  CD_CD4T_degs = CD_CD4T_degs,
  UC_CD8T_degs = UC_CD8T_degs,
  CD_CD8T_degs = CD_CD8T_degs,
  UC_ILC_degs = UC_ILC_degs,
  CD_ILC_degs = CD_ILC_degs,
  UC_B_Plasma_degs = UC_B_Plasma_degs,
  CD_B_Plasma_degs = CD_B_Plasma_degs
)

write_rds(dataframe_map, "all_degs.rds")

dataframe_map <- lapply(dataframe_map, function(x) {
  x <- x %>% filter(p_val_adj < 0.05)
  return(x)
})

all_genes = [[...]]
all_genes_results <- list()

for (name in names(dataframe_map)) {
  all_genes_results[[paste0(name, "_all_genes")]] <- data.frame(gene = character(), trait = character(), cell_type = character(), cluster = character(), stringsAsFactors = FALSE)
}

for (name in names(dataframe_map)) {
  all_genes_results[[paste0(name, "_all_genes")]] <- data.frame(
    gene = character(), trait = character(), cell_type = character(), cluster = character(),
    p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(),
    stringsAsFactors = FALSE
  )
}

for (i in 1:nrow(all_genes)) {
  gene <- all_genes$gene[i]
  trait <- all_genes$trait[i]
  cluster <- all_genes$cluster[i]
  
  if (trait == "IBD") {
    # if trait is IBD, search in both UC and CD for the given cluster
    uc_dataframe_name <- paste("UC", cluster, "degs", sep = "_")
    cd_dataframe_name <- paste("CD", cluster, "degs", sep = "_")
    
    # check and add gene to UC dataframe if it exists
    if (uc_dataframe_name %in% names(dataframe_map)) {
      uc_df <- dataframe_map[[uc_dataframe_name]]
      if (gene %in% rownames(uc_df)) {
        gene_data <- uc_df[gene, ]
        all_genes_results[[paste0(uc_dataframe_name, "_all_genes")]] <- rbind(
          all_genes_results[[paste0(uc_dataframe_name, "_all_genes")]],
          data.frame(
            gene = gene, trait = "UC", cell_type = all_genes$cell_type[i], cluster = cluster,
            p_val = gene_data$p_val, avg_log2FC = gene_data$avg_log2FC, pct.1 = gene_data$pct.1, pct.2 = gene_data$pct.2, p_val_adj = gene_data$p_val_adj,
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    # check and add gene to CD dataframe if it exists
    if (cd_dataframe_name %in% names(dataframe_map)) {
      cd_df <- dataframe_map[[cd_dataframe_name]]
      if (gene %in% rownames(cd_df)) {
        gene_data <- cd_df[gene, ]
        all_genes_results[[paste0(cd_dataframe_name, "_all_genes")]] <- rbind(
          all_genes_results[[paste0(cd_dataframe_name, "_all_genes")]],
          data.frame(
            gene = gene, trait = "CD", cell_type = all_genes$cell_type[i], cluster = cluster,
            p_val = gene_data$p_val, avg_log2FC = gene_data$avg_log2FC, pct.1 = gene_data$pct.1, pct.2 = gene_data$pct.2, p_val_adj = gene_data$p_val_adj,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  } else {
    # for non-IBD traits, use the original code logic
    dataframe_name <- paste(trait, cluster, "degs", sep = "_")
    
    # check if the dataframe exists in the map
    if (dataframe_name %in% names(dataframe_map)) {
      df <- dataframe_map[[dataframe_name]]
      
      # check if the gene is in the rownames of the dataframe
      if (gene %in% rownames(df)) {
        gene_data <- df[gene, ]
        # append the gene to the corresponding all_genes dataframe with expression data
        all_genes_results[[paste0(dataframe_name, "_all_genes")]] <- rbind(
          all_genes_results[[paste0(dataframe_name, "_all_genes")]],
          data.frame(
            gene = gene, trait = trait, cell_type = all_genes$cell_type[i], cluster = cluster,
            p_val = gene_data$p_val, avg_log2FC = gene_data$avg_log2FC, pct.1 = gene_data$pct.1, pct.2 = gene_data$pct.2, p_val_adj = gene_data$p_val_adj,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

verified_genes <- do.call(rbind, all_genes_results)

verified_genes <- verified_genes %>% 
  tibble::rownames_to_column(var = "data_id") %>%
  select(-data_id, -cell_type, -pct.1, -pct.2, p_val) %>%
  distinct()

print(length(unique(verified_genes$gene)))

write.xlsx(verified_genes, "..._genes_verified.xlsx")
write_rds(all_genes_results, "..._genes_verified.rds")
