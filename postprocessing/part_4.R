setwd("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_4/")

library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(pheatmap)

scIBD_seurat <- readRDS("scIBD_filtered.rds")
scIBD_degs <- readRDS("all_degs.rds")

# bar plot for number of cells in each cluster in each sample
cell_counts_df <- as.data.frame(cell_counts)

total_cells <- sum(cell_counts_df$Cell_Count)

cell_counts <- table(scIBD_seurat@meta.data$sample, scIBD_seurat@meta.data$major_cluster)

colnames(cell_counts_df) <- c("Sample", "Major_Cluster", "Cell_Count")

cell_counts_df <- cell_counts_df %>%
  filter(Cell_Count > 0)

ggplot(cell_counts_df, aes(x = Sample, y = Cell_Count, fill = Major_Cluster)) +
  geom_bar(stat = "identity") +  
  labs(title = "Number of Cells in Each Major Cluster by Sample",
       x = "Sample",
       y = "Number of Cells",
       fill = "Major Cluster") +  
  scale_fill_brewer(palette = "Set3") +  
  theme_minimal() +               
  theme(
    text = element_text(color = "black"),
    
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(color = "black"),
    
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    
    plot.title = element_text(color = "black", hjust = 0.5)
  )

ggsave("cell_counts_sample.png", width = 10, height = 6, dpi = 300)
write.csv(cell_counts_df, "cell_counts_sample.csv", row.names = FALSE)

# calculating the total number of cells in each major cluster in total
total_cells <- sum(cell_counts_df$Cell_Count)

cell_counts_df <- cell_counts_df %>%
  group_by(Major_Cluster) %>%
  summarise(Total_Cells = sum(Cell_Count))

# bar plot for number of cells in each cluster in each condition
conditions <- c("UC_PBMC", "Healthy_PBMC", "CD_PBMC")
filtered_data <- scIBD_seurat@meta.data %>%
  filter(disease %in% conditions)

cell_counts <- table(filtered_data$disease, filtered_data$major_cluster)

cell_counts_df <- as.data.frame(cell_counts)

colnames(cell_counts_df) <- c("Condition", "Major_Cluster", "Cell_Count")

cell_counts_df <- cell_counts_df %>%
  filter(Cell_Count > 0)

ggplot(cell_counts_df, aes(x = Condition, y = Cell_Count, fill = Major_Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Number of Cells in Each Major Cluster by Condition",
       x = "Condition",
       y = "Number of Cells",
       fill = "Major Clusters") +
  scale_fill_brewer(palette = "Set3") + 
  theme_minimal() +
  theme(
    text = element_text(color = "black"),
    
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"), 
    axis.text.y = element_text(color = "black"),

    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),

    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    
    plot.title = element_text(color = "black", hjust = 0.5)
  )

ggsave("cell_counts_conditions.png", width = 10, height = 6, dpi = 300)
write.csv(cell_counts_df, "cell_counts_conditions.csv", row.names = FALSE)

# UMAP
umap_data <- scIBD_seurat@meta.data

ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = major_cluster)) +
  geom_point(size = 1, alpha = 0.8) +  
  theme_minimal() +  
  labs(title = "UMAP Plot Colored by Major Cluster",
       x = "UMAP_1", y = "UMAP_2", color = "Major Clusters") +  
  theme(
    text = element_text(color = "black"),
    
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) +
  scale_color_brewer(palette = "Set3")

ggsave("umap_major_clusters_manual.png", width = 8, height = 6, dpi = 300)

# the number of DEGs each major cluster 
scIBD_degs <- lapply(scIBD_degs, function(x) {
  x <- x %>% filter(p_val_adj < 0.05)
  return(x)
})

uc_genes <- unique(c(rownames(scIBD_degs$UC_Myeloid_degs), rownames(scIBD_degs$UC_CD4T_degs),
                  rownames(scIBD_degs$UC_CD8T_degs), rownames(scIBD_degs$UC_ILC_degs),
                  rownames(scIBD_degs$UC_B_Plasma_degs)))

cd_genes <- unique(c(rownames(scIBD_degs$CD_Myeloid_degs), rownames(scIBD_degs$CD_CD4T_degs),
                  rownames(scIBD_degs$CD_CD8T_degs), rownames(scIBD_degs$CD_ILC_degs),
                  rownames(scIBD_degs$CD_B_Plasma_degs)))

genes <- unique(c(uc_genes, cd_genes))

length(intersect(uc_genes, cd_genes))

# heatmap for verified prioritized genes
scIBD_ids <- read_excel("scIBD_ids.xlsx", sheet = "Sheet2")
scIBD_seurat <- readRDS("scIBD_filtered.rds")
genes_verified <- read_excel("prioritized_genes_verified.xlsx")

genes_of_interest <- unique(genes_verified$gene)
seurat_subset <- subset(scIBD_seurat, features = genes_of_interest)

expression_matrix <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")

metadata <- seurat_subset@meta.data

combined_data <- data.frame(metadata, t(as.matrix(expression_matrix)))

avg_expression <- combined_data %>%
  select(-c(sample, subject, study, nUMI, nGene, rp_pct, mt_pct, minor_cluster, UMAP_1, UMAP_2, TSNE_1, TSNE_2, gUMAP_1
            , gUMAP_2, gTSNE_1, gTSNE_2, stage, tissue, tissue.sub, nCount_RNA, nFeature_RNA))

avg_expression <- avg_expression %>%
  group_by(major_cluster, disease) %>%
  summarize(across(3:41, mean))  

heatmap_matrix <- as.matrix(avg_expression[,-c(1,2)])

rownames(heatmap_matrix) <- paste(avg_expression$major_cluster, avg_expression$disease, sep = "_")

tr_heatmap_matrix <- t(heatmap_matrix)

ann_df <- data.frame(
  major_cluster = avg_expression$major_cluster,
  disease = avg_expression$disease
)

ann_df$disease <- sapply(ann_df$disease, function(x) { strsplit(as.character(x), "_")[[1]][1]})
ann_df$disease <- factor(ann_df$disease)

rownames(ann_df) <- colnames(tr_heatmap_matrix)

annotation_colors <- list(
  major_cluster = c(
    "Myeloid" = "#9CD1C7",
    "CD4T" = "#FFFFB4",
    "CD8T" = "#BDBAD7",
    "ILC" = "#EB8677",
    "B_Plasma" = "#8AB0D0"
  ),
  disease = c(
    "UC" = "#E93423",
    "CD" = "#306FBA",
    "Healthy" = "grey"
  )
)

png("heatmap_of_pgenes_complex.png", width = 10, height = 10, units = "in", res = 300)

pheatmap(
  tr_heatmap_matrix,
  cluster_rows = TRUE,   # Cluster genes
  cluster_cols = FALSE,  # Do not cluster columns
  show_rownames = TRUE,  # Show gene names (rows)
  show_colnames = FALSE,  # Show cluster/condition names (columns)
  fontsize_col = 10,     # Adjust column font size for readability
  fontsize_row = 8,      # Adjust row font size for readability
  color = colorRampPalette(c("white", "brown"))(100),
  annotation_colors = annotation_colors,  # Custom color scale
  annotation_col = ann_df,  # Add disease annotation
  name = "Expression Level"  # Set the legend title
)

dev.off()











