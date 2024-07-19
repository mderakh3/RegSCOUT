setwd("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_2")

library(stringr)
library(tidyverse)
library(tidyr)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)

# preapring table 1 of the paper
risk_regions_ppa <- read_excel("risk_regions_ppa.xlsx")

risk_regions_ppa <- risk_regions_ppa %>% 
  separate_rows(cell, sep = ",") %>%
  distinct()

table_1 <- data.frame(matrix(ncol = 5, nrow = length(unique(risk_regions_ppa$cell))))
colnames(table_1) <- c("cell", "total_rmp",  "unique_rmp", "total_po_snp", "total_tf")

table_1$cell <- unique(risk_regions_ppa$cell)

for (i in 1:nrow(table_1)) {
  table_1$total_rmp[i] <- length(unique(risk_regions_ppa$region[risk_regions_ppa$cell == table_1$cell[i]]))
}

specific_regions <- risk_regions_ppa %>%
  group_by(region) %>%
  filter(n() == 1) %>%
  ungroup()

cell_specific_counts <- specific_regions %>%
  group_by(cell) %>%
  summarize(count = n())

for (i in 1:nrow(table_1)) {
  if (table_1$cell[i] %in% cell_specific_counts$cell) {
    table_1$unique_rmp[i] <- cell_specific_counts$count[cell_specific_counts$cell == table_1$cell[i]]
  } else {
    table_1$unique_rmp[i] <- 0
  }
}

rm(specific_regions, cell_specific_counts)

risk_regions_ratio <- read_excel("risk_regions_ratio.xlsx")

risk_regions_ratio$snp = gsub(".*-", "", risk_regions_ratio$TFSNP)
risk_regions_ratio$tf <- risk_regions_ratio$TFSNP

for (i in 1:nrow(risk_regions_ratio)) {
  current_string <- risk_regions_ratio$TFSNP[i]
  dash_count <- sum(strsplit(current_string, "")[[1]] == "-")
  if (dash_count >= 2) {
    risk_regions_ratio$tf[i] <- gsub("^(.*?-.*?)-.*$", "\\1", current_string)
  } else if (dash_count == 1) {
    risk_regions_ratio$tf[i] <- gsub("^(.*?)-.*$", "\\1", current_string)
  }
}

risk_regions_ratio <- risk_regions_ratio %>% 
  separate_rows(cell, sep = ",") %>%
  distinct()

for (i in 1:nrow(table_1)) {
  cell <- table_1$cell[i]
  table_1$total_po_snp[i] <- length(unique(risk_regions_ratio$snp[risk_regions_ratio$cell == cell]))
  table_1$total_tf[i] <- length(unique(risk_regions_ratio$tf[risk_regions_ratio$cell == cell]))
}

write.table(table_1, "table_1.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# heatmap for rmp promoters and genes in cell types
rmp_promoter_data <- read_excel("rmp_promoter_data.xlsx")

rmp_promoter_data <- rmp_promoter_data %>% 
  dplyr::select(rmp_prom, gene_name) %>%
  distinct()

gene_prom_matrix <- merge(rmp_promoter_data, risk_regions_ppa, by.x = "rmp_prom", by.y = "region", all = T) %>%
  dplyr::select(-SNP, -sumPPA) %>%
  filter(!is.na(gene_name)) %>%
  dplyr::select(-rmp_prom) %>%
  distinct() %>%
  mutate(value = 1) %>%
  spread(cell, value, fill = 0) %>%
  column_to_rownames(var = "gene_name")

ht <- Heatmap(
  as.matrix(gene_prom_matrix), 
  name = "Gene Presence",
  col = c("0" = "white", "1" = "red"),
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(2, "cm"),
  border = TRUE, # Add borders around the entire heatmap
  rect_gp = gpar(col = "black", lwd = 1), # Add borders around the tiles
  row_title = "Genes",
  column_title = "Immune Cell Subtypes",
  row_title_gp = gpar(fontsize = 15), # Increase font size of the y-axis title
  column_title_gp = gpar(fontsize = 15) # Increase font size of the x-axis title
)

png(filename = "prom_gene_heatmap.png", width = 1800, height = 2400, res = 300)
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# heatmap for rmp enhancers and genes in cell types
cic_rmp_enhancers <- read_excel("cic_rmp_enhancers.xlsx")

cic_rmp_enhancers <- cic_rmp_enhancers %>% select(-Peak1, -Peak2, -type, -prom_side, -promoter)
cic_rmp_enhancers[is.na(cic_rmp_enhancers)] <- 0

cic_rmp_enhancers <- cic_rmp_enhancers %>% mutate(b = ifelse(b > 0, 1, 0), t = ifelse(t > 0, 1, 0),
                                                  nk = ifelse(nk > 0, 1, 0), mono = ifelse(mono > 0, 1, 0))
cic_rmp_enhancers <- cic_rmp_enhancers %>% select(-enhancer)
cic_rmp_enhancers <- cic_rmp_enhancers %>% separate_rows(genes, sep = ",")
cic_rmp_enhancers <- cic_rmp_enhancers %>% select(genes, b, t, nk, mono)

gene_enh_matrix <- data.frame(matrix(nrow = length(unique(cic_rmp_enhancers$genes)), ncol = ncol(cic_rmp_enhancers) - 1))
rownames(gene_enh_matrix) <- unique(cic_rmp_enhancers$genes)
colnames(gene_enh_matrix) <- colnames(cic_rmp_enhancers)[-1]

for (i in 1:nrow(gene_enh_matrix)) {
  for (j in 1:ncol(gene_enh_matrix)) {
    gene_enh_matrix[i, j] <- ifelse(any(cic_rmp_enhancers$genes == rownames(gene_enh_matrix)[i] & cic_rmp_enhancers[[j + 1]] == 1), 1, 0)
  }
}

colnames(gene_enh_matrix) <- c("B_cells", "T_cells", "NK_cells", "Monocytes")

ht <- Heatmap(
  as.matrix(gene_enh_matrix), 
  name = "Gene Presence",
  col = c("0" = "white", "1" = "red"),
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 1), # Increase font size of row names
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(2, "cm"),
  border = TRUE, # Add borders around the entire heatmap
  rect_gp = gpar(col = "black", lwd = 0.1), # Add borders around the tiles
  row_title = "Genes",
  column_title = "Immune Cell Types",
  row_title_gp = gpar(fontsize = 25), # Increase font size of the y-axis title
  column_title_gp = gpar(fontsize = 25), # Increase font size of the x-axis title
  column_names_gp = gpar(fontsize = 15), # Increase font size of the column names
  heatmap_legend_param = list(
    legend_gp = gpar(fontsize = 25) # Increase font size of legend
  )
)

png(filename = "enh_gene_heatmap.png", width = 1600, height = 3600, res = 300)
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# venn diagram for shared rmps 
cic_rmp_enhancers <- read_excel("cic_rmp_enhancers.xlsx")

rmp_enhancers <- unique(cic_rmp_enhancers$enhancer)
rmp_promoters <- unique(rmp_promoter_data$rmp_prom)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(rmp_enhancers, rmp_promoters),
  category.names = c("Potential Enhancer RMPs", "Potential Promoter RMPs"),
  filename = "venn_diagram.png",
  output = TRUE,
  imagetype = "png",
  main = "Shared RMPs"
)

png(filename = "venn_diagram.png", width = 1200, height = 1200, res = 300)
dev.off()

# hic verfied loops table
cic_rmp_enhancers <- read_excel("cic_rmp_enhancers.xlsx")
cic_rmp_enhancers <- cic_rmp_enhancers %>% select(-Peak1, -Peak2, -type, -prom_side)

cic_data_names <- c("b", "t", "nk", "mono")
cell_cic_data <- list()

for (i in 1:length(cic_data_names)) {
  cell_cic_data[[i]] <- cic_rmp_enhancers %>% select(promoter, enhancer, genes, cic_data_names[i])
  cell_cic_data[[i]] <- cell_cic_data[[i]] %>% filter(!is.na(cell_cic_data[[i]][[4]]))
  cell_cic_data[[i]] <- cell_cic_data[[i]] %>% select(-cic_data_names[i])
}

names(cell_cic_data) <- cic_data_names

hic_data_names <- c("b_pcHiC_Cic", "mono_pcHiC_Cic", "nk_intHiC_Cic", "t_pcHiC_Cic")
cell_hic_data <- list()

for (i in 1:length(hic_data_names)) {
  cell_hic_data[[i]] <- read.table(paste0(hic_data_names[i], ".txt"), header = T, sep = "\t")
}

names(cell_hic_data) <- hic_data_names

for (i in 1:length(cell_hic_data)) {
  cell_hic_data[[i]] <- cell_hic_data[[i]] %>% select(cic_promoter, cic_enhancer, cic_genes)
}

for (i in 1:length(cell_cic_data)) {
  cell_cic_data[[i]]$loop <- paste(cell_cic_data[[i]]$promoter, cell_cic_data[[i]]$enhancer, sep = "_")
  cell_cic_data[[i]] <- cell_cic_data[[i]] %>% select(-promoter, -enhancer)
  cell_cic_data[[i]] <- cell_cic_data[[i]] %>% distinct()
}

for (i in 1:length(cell_hic_data)) {
  cell_hic_data[[i]]$loop <- paste(cell_hic_data[[i]]$cic_promoter, cell_hic_data[[i]]$cic_enhancer, sep = "_")
  cell_hic_data[[i]] <- cell_hic_data[[i]] %>% select(-cic_promoter, -cic_enhancer)
  cell_hic_data[[i]] <- cell_hic_data[[i]] %>% distinct()
}

loop_data <- data.frame(cell_type = character(), cicero_loops = numeric(), 
                        hic_loops = numeric(), shared_percentage = numeric())

for (i in 1:length(cell_cic_data)) {
  cicero_loops <- nrow(cell_cic_data[[i]])
  hic_loops <- nrow(cell_hic_data[[i]])
  verification_percentage <- hic_loops / cicero_loops * 100
  
  loop_data <- rbind(loop_data, data.frame(cell_type = cic_data_names[i], cicero_loops = cicero_loops, 
                                           hic_loops = hic_loops, 
                                           verification_percentage = verification_percentage))
}

loop_data <- loop_data %>% arrange(desc(verification_percentage))

write.table(loop_data, "loop_data.txt", sep = "\t", row.names = F)

# cicero and hic heatmap
cic_rmp_enhancers <- read_excel("cic_rmp_enhancers.xlsx")
cic_rmp_enhancers <- cic_rmp_enhancers %>% select(-Peak1, -Peak2, -type, -prom_side)

hic_data_names <- c("b_pcHiC_Cic", "mono_pcHiC_Cic", "nk_intHiC_Cic", "t_pcHiC_Cic")
cell_hic_data <- list()

for (i in 1:length(hic_data_names)) {
  cell_hic_data[[i]] <- read.table(paste0(hic_data_names[i], ".txt"), header = TRUE, sep = "\t")
}

names(cell_hic_data) <- hic_data_names

cicero_heatmap <- function(cic_rmp_enhancers, cell_type, hic_data) {
  cell_cic_data <- cic_rmp_enhancers %>% select(promoter, enhancer, cell_type)
  cell_cic_data <- cell_cic_data[!is.na(cell_cic_data[[cell_type]]), ]
  cell_cic_data <- cell_cic_data %>% distinct(promoter, enhancer, .keep_all = TRUE)
  cell_cic_matrix <- cell_cic_data %>% pivot_wider(names_from = enhancer, values_from = cell_type)
  cell_cic_matrix[is.na(cell_cic_matrix)] <- 0
  cell_cic_matrix <- tibble::column_to_rownames(cell_cic_matrix, "promoter")
  cell_cic_matrix <- cell_cic_matrix[rowSums(cell_cic_matrix) != 0, ]
  cell_cic_matrix[cell_cic_matrix > 0] <- 1
  
  confirmed <- as.matrix(cell_cic_matrix)   # Mark Hi-C confirmed regions
  for (i in 1:nrow(hic_data)) {
    prom <- hic_data$cic_promoter[i]
    enh <- hic_data$cic_enhancer[i]
    if (prom %in% rownames(confirmed) && enh %in% colnames(confirmed)) {
      confirmed[prom, enh] <- 2
    }
  }
  return(list(matrix = cell_cic_matrix, confirmed = confirmed))
}

cell_types <- c("b", "mono", "nk", "t")
hic_data_map <- list(
  b = "b_pcHiC_Cic",
  mono = "mono_pcHiC_Cic",
  nk = "nk_intHiC_Cic",
  t = "t_pcHiC_Cic"
)
cell_cic_matrices <- list()

for (i in 1:length(cell_types)) {
  hic_data <- cell_hic_data[[hic_data_map[[cell_types[i]]]]]
  cell_cic_matrices[[cell_types[i]]] <- cicero_heatmap(cic_rmp_enhancers, cell_types[i], hic_data)
}

heatmap_function <- function(cell_cic_matrices, cell_type) {
  cell_type_names <- list(
    b = "B",
    t = "T",
    nk = "NK",
    mono = "Mono"
  )
  
  cell_type_full <- cell_type_names[[cell_type]]
  cell_cic_matrix <- cell_cic_matrices[[cell_type]]$matrix
  confirmed <- cell_cic_matrices[[cell_type]]$confirmed
  
  ht <- Heatmap(
    confirmed,
    name = "Regions",
    col = c("0" = "grey", "1" = "red", "2" = "blue"),
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 3),
    row_names_side = "right",
    row_dend_width = unit(3, "cm"),
    column_dend_height = unit(2, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "black", lwd = 1),
    row_title = "Promoter-Genes",
    row_title_gp = gpar(fontsize = 45),
    column_title = "Distal RMPs",
    column_title_gp = gpar(fontsize = 45),
    heatmap_legend_param = list(
      at = c(0, 1, 2),
      labels = c("Not Co-accessible", "Co-accessible", "Hi-C Verified"),
      title_gp = gpar(fontsize = 35, fontface = "bold"),  # Increase title font size
      labels_gp = gpar(fontsize = 25)  # Increase labels font size
    )
  )
  
  png(filename = paste0(cell_type, "_cic_heatmap.png"), width = 7200, height = 4800, res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grid.text(cell_type_full, x = unit(0.04, "npc"), y = unit(0.97, "npc"), gp = gpar(fontsize = 40, fontface = "bold"))
  dev.off()
}

for (i in 1:length(cell_types)) {
  heatmap_function(cell_cic_matrices, cell_types[i])
}

# hic confirmed genes heatmap
hic_data_names <- c("b_pcHiC_Cic", "mono_pcHiC_Cic", "nk_intHiC_Cic", "t_pcHiC_Cic")
cell_hic_data <- list()

for (i in 1:length(hic_data_names)) {
  cell_hic_data[[i]] <- read.table(paste0(hic_data_names[i], ".txt"), header = TRUE, sep = "\t")
}

names(cell_hic_data) <- hic_data_names

separate_genes <- function(data, cic_col, pchic_col, cic_sep = ",", pchic_sep = ";") {
  data <- data %>%
    mutate(cic_genes = strsplit(as.character(!!sym(cic_col)), cic_sep)) %>%
    unnest(cic_genes) %>%
    mutate(pchic_genes = strsplit(as.character(!!sym(pchic_col)), pchic_sep)) %>%
    unnest(pchic_genes) %>%
    filter(cic_genes == pchic_genes)
  return(data)
}

cell_hic_data$t_pcHiC_Cic <- separate_genes(cell_hic_data$t_pcHiC_Cic, "cic_genes", "pchic_ba_genes")
cell_hic_data$b_pcHiC_Cic <- separate_genes(cell_hic_data$b_pcHiC_Cic, "cic_genes", "pchic_ba_genes")
cell_hic_data$mono_pcHiC_Cic <- separate_genes(cell_hic_data$mono_pcHiC_Cic, "cic_genes", "pchic_ba_genes")

separate_nk_genes <- function(data, cic_col, cic_sep = ",") {
  data <- data %>%
    mutate(cic_genes = strsplit(as.character(!!sym(cic_col)), cic_sep)) %>%
    unnest(cic_genes)
  return(data)
}

cell_hic_data$nk_intHiC_Cic <- separate_nk_genes(cell_hic_data$nk_intHiC_Cic, "cic_genes")

get_unique_genes <- function(data, gene_col) { # Get unique genes for each cell type
  unique(data[[gene_col]])
}

genes_t <- get_unique_genes(cell_hic_data$t_pcHiC_Cic, "cic_genes")
genes_b <- get_unique_genes(cell_hic_data$b_pcHiC_Cic, "cic_genes")
genes_mono <- get_unique_genes(cell_hic_data$mono_pcHiC_Cic, "cic_genes")
genes_nk <- get_unique_genes(cell_hic_data$nk_intHiC_Cic, "cic_genes")

all_genes <- unique(c(genes_t, genes_b, genes_mono, genes_nk))

gene_matrix <- data.frame( # Initialize a binary matrix
  Gene = all_genes,
  T_cells = 0,
  B_cells = 0,
  Monocytes = 0,
  NK_Cells = 0,
  stringsAsFactors = FALSE
)

update_binary_matrix <- function(gene_matrix, genes, cell_type_col) {
  gene_matrix[[cell_type_col]] <- ifelse(gene_matrix$Gene %in% genes, 1, 0)
  return(gene_matrix)
}

gene_matrix <- update_binary_matrix(gene_matrix, genes_t, "T_cells")
gene_matrix <- update_binary_matrix(gene_matrix, genes_b, "B_cells")
gene_matrix <- update_binary_matrix(gene_matrix, genes_mono, "Monocytes")
gene_matrix <- update_binary_matrix(gene_matrix, genes_nk, "NK_Cells")
gene_matrix <- tibble::column_to_rownames(gene_matrix, var = "Gene")

ht <- Heatmap(
  as.matrix(gene_matrix), 
  name = "Gene Presence",
  col = c("0" = "white", "1" = "red"),
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(2, "cm"),
  border = TRUE, # Add borders around the entire heatmap
  rect_gp = gpar(col = "black", lwd = 1), # Add borders around the tiles
  row_title = "Genes",
  column_title = "Immune Cell Types",
  row_title_gp = gpar(fontsize = 15), # Increase font size of the y-axis title
  column_title_gp = gpar(fontsize = 15) # Increase font size of the x-axis title
)

png(filename = "hic_gene_heatmap.png", width = 1800, height = 2400, res = 300)
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# eqtl table
files <- list.files(path = "/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_2/eQTL", 
                    pattern = "*.txt", full.names = TRUE)

process_file <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$external_gene_name[df$external_gene_name == ""] <- df$gene_id[df$external_gene_name == ""]
  df_combined <- df %>%
    group_by(rsid) %>%
    summarise(external_gene_name = paste(unique(external_gene_name), collapse = ", "),
              gene_id = paste(unique(gene_id), collapse = ", "),
              pvalue = min(pvalue)) %>% # Optional: take the minimum p-value
    ungroup()
  return(df_combined)
}

eQTL_data <- lapply(files, process_file)
eQTL_combined <- bind_rows(eQTL_data)
names(eQTL_data) <- basename(files)

length(unique(eQTL_combined$rsid))
genes <- unique(eQTL_combined %>% separate_rows(external_gene_name, sep = ", ") %>% pull(external_gene_name))

complete_table <- read.csv("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_2/complete_table.csv")
complete_table <- complete_table[!is.na(complete_table$eqtl_gene),]
complete_table <- complete_table %>% separate_rows(eqtl_gene, sep = ", ")
eqtl_rmp <- unique(complete_table$rmp)

table_3 <- data.frame(matrix(ncol = 4, nrow = 13))
colnames(table_3) <- c("cell_type", "eqtl_rmp", "eqtl_snps", "eqtl_genes")
cell_types <- c("cMono", "ncMono", "iMono", "adaptive_NK", "cyto_nk", 
                "tReg", "mem_b", "naive_b", "act_cd4_t", "naive_cd4_t", "naive_cd8_t", 
                "mem_cd8_t", "cyto_cd8_t")
table_3$cell_type <- cell_types

for (i in 1:nrow(table_3)) {
  cell_type <- table_3$cell_type[i]
  com_table <- complete_table[complete_table$cell_type == cell_type,]
  table_3$eqtl_rmp[i] <- length(unique(com_table$rmp))
  table_3$eqtl_snps[i] <- length(unique(com_table$effect_snp))
  table_3$eqtl_genes[i] <- length(unique(com_table$eqtl_gene))
}

write.table(table_3, "supplementary_table_2c.txt", row.names = FALSE, sep = "\t")

# eqtl heatmap
complete_table <- complete_table %>% select(eqtl_gene, cell_type)
complete_table <- complete_table %>% distinct()

all_genes <- unique(complete_table$eqtl_gene)

eqtl_matrix <- data.frame( # Initialize a binary matrix
  Gene = all_genes,
  cMono = 0,
  ncMono = 0,
  iMono = 0,
  adaptive_NK = 0,
  cyto_nk = 0,
  tReg = 0,
  mem_b = 0,
  naive_b = 0,
  act_cd4_t = 0,
  naive_cd4_t = 0,
  naive_cd8_t = 0,
  mem_cd8_t = 0,
  cyto_cd8_t = 0,
  stringsAsFactors = FALSE
)

eqtl_matrix <- tibble::column_to_rownames(eqtl_matrix, var = "Gene")

for (i in 1:nrow(complete_table)) {
  gene <- complete_table$eqtl_gene[i]
  cell_type <- complete_table$cell_type[i]
  eqtl_matrix[gene, cell_type] <- 1
}

ht <- Heatmap(
  as.matrix(eqtl_matrix), 
  name = "Gene Presence",
  col = c("0" = "white", "1" = "red"),
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 3),
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(2, "cm"),
  border = TRUE, # Add borders around the entire heatmap
  rect_gp = gpar(col = "black", lwd = 0.1), # Add borders around the tiles
  row_title = "Genes",
  column_title = "Immune Cell Subtypes",
  row_title_gp = gpar(fontsize = 25), # Increase font size of the y-axis title
  column_title_gp = gpar(fontsize = 25), # Increase font size of the x-axis title
  column_names_gp = gpar(fontsize = 15), # Increase font size of the column names
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 15)  # Increase labels font size
  )
)

png(filename = "eqtl_gene_heatmap.png", width = 1800, height = 3600, res = 300)
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

# all 
complete_table <- read.csv("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_2/complete_table.csv")

complete_table[complete_table == "not found"] <- NA
complete_table[complete_table == "not applicable"] <- NA

geneless_rows <- complete_table %>% filter(is.na(directly_mapped_gene) & is.na(cicero_gene) & is.na(hic_gene) & is.na(eqtl_gene))
length(unique(geneless_rows$rmp))

complete_table <- complete_table %>% filter(!is.na(directly_mapped_gene) | !is.na(cicero_gene) | !is.na(hic_gene) | !is.na(eqtl_gene))
length(unique(complete_table$rmp))

# disgennet
complete_table <- read.csv("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_2/complete_table.csv")

complete_table$gene <- NA

for (i in 1:nrow(complete_table)) {
  complete_table$gene[i] <- paste0(complete_table$directly_mapped_gene[i], ", ", complete_table$cicero_gene[i], ", ", complete_table$hic_gene[i], ", ", complete_table$eqtl_gene[i])
}

complete_table <- complete_table %>% select(-directly_mapped_gene, -cicero_gene, -hic_gene, -eqtl_gene)
complete_table <- complete_table %>% separate_rows(gene, sep = ", ")
complete_table$gene[complete_table$gene == "NA"] <- NA
complete_table$gene[complete_table$gene == "not found"] <- NA
complete_table <- complete_table[!is.na(complete_table$gene),]

complete_genes <- complete_table %>% pull(gene)
complete_genes <- unique(complete_genes)

complete_table <- complete_table %>% separate_rows(tf, sep = ", ")

prioritized_table <- read.table("prioritized_table.txt" , header = TRUE, sep = "\t")
prioritized_table <- prioritized_table %>% separate_rows(gene_score, sep = ", ")

prioritized_table$gene_score <- gsub("\\(.*\\)", "", prioritized_table$gene_score)
prioritized_table$gene_score <- gsub(" ", "", prioritized_table$gene_score)

prioritized_genes <- unique(prioritized_table$gene_score)

library(devtools) # retrieving information from DisGenNet database
library(disgenet2r)

api_key <- "9419340c-78da-4d4f-9020-3bab3304d6e1"
Sys.setenv(DISGENET_API_KEY= api_key)


results <- gene2disease(gene = prioritized_genes, vocabulary = "HGNC",
                         database = "ALL", score = c(0.1,1))

write_rds(results, "86_genes_to_diseases.rds")

diseasesOfInterest <- paste0("UMLS_",c("C0021390", "C0010346", "C0009324"))

all_disgennet_genes <- disease2gene(
  disease = diseasesOfInterest,
  database = "ALL",
  score =c(0.1,1),
  verbose  = TRUE )

write_rds(all_disgennet_genes, "all_disgennet_genes.rds")

curated_disgennet_genes <- disease2gene(
  disease = diseasesOfInterest,
  database = "CURATED",
  score =c(0.1,1),
  verbose  = TRUE )

disgennet_df <- unique(curated_disgennet_genes@qresult[  ,c("gene_symbol", "disease_name","score", "chemicalsIncludedInEvidence", "yearFinal", "diseaseUMLSCUI")] )

disgennet_df$chemicals <- NA
for (i in 1:nrow(disgennet_df)) {
  if (length(disgennet_df$chemicalsIncludedInEvidence[[i]]) > 0) {
    disgennet_df$chemicals[i] <- paste(disgennet_df$chemicalsIncludedInEvidence[[i]], collapse = ", ")
  }
}

write.table(disgennet_df, "disgennet_df.txt", row.names = FALSE, sep = "\t")

all_disgennet_genes <- readRDS("all_disgennet_genes.rds")
curated_disgennet_genes <- readRDS("curated_disgennet_genes.rds")

all_disgennet_genes <- unique(all_disgennet_genes@qresult$gene_symbol)
curated_disgennet_genes <- unique(curated_disgennet_genes@qresult$gene_symbol)

gene_sets <- list(
  All_DisGeNET_Genes = all_disgennet_genes,
  Curated_DisGeNET_Genes = curated_disgennet_genes,
  Prioritized_Genes = prioritized_genes,
  Complete_Genes = complete_genes
)

gene_sets_named <- list(
  "All DisGeNET Genes" = all_disgennet_genes,
  "Curated DisGeNET Genes" = curated_disgennet_genes,
  "Prioritized Genes" = prioritized_genes,
  "Complete Genes" = complete_genes
)

library(eulerr)
fit <- euler(gene_sets_named)

venn.plot = plot(fit, 
     fills = list(fill = c("cornflowerblue", "green", "yellow", "darkorchid1"), alpha = 0.5), 
     edges = TRUE,
     quantities = TRUE,
     labels = list(cex = 1.5),
     legend = list(labels = names(gene_sets_named), side = "right", cex = 1.5))

png(filename = "venn_diagram.png", width = 4800, height = 4800, res = 300)
grid.draw(venn.plot)
dev.off()

prioritized_genes <- setdiff(prioritized_genes, c(all_disgennet_genes, curated_disgennet_genes))

# prioritized genes, tfs, snps, and loci
prioritized_genes <- unique(prioritized_table$gene_score)
prioritized_cells <- unique(prioritized_table$cell_type)
gene_matrix <- matrix(0, nrow = length(prioritized_genes), ncol = length(prioritized_cells))
rownames(gene_matrix) <- prioritized_genes
colnames(gene_matrix) <- prioritized_cells

for (i in 1:nrow(prioritized_table)) { # gene_matrix
  gene_matrix[prioritized_table$gene_score[i], prioritized_table$cell_type[i]] <- 1
}

prioritized_tfs_table <- read.table("prioritized_table.txt", header = TRUE, sep = "\t")
prioritized_tfs_table <- prioritized_tfs_table %>% separate_rows(tf, sep = ", ")
prioritized_tfs_table <- distinct(prioritized_tfs_table)
prioritized_tfs_table$tf <- gsub("\\(.*\\)", "", prioritized_tfs_table$tf)
prioritized_tfs <- unique(prioritized_tfs_table$tf)

tf_matrix <- matrix(0, nrow = length(prioritized_tfs), ncol = length(prioritized_cells))
rownames(tf_matrix) <- prioritized_tfs
colnames(tf_matrix) <- prioritized_cells


for (i in 1:nrow(prioritized_tfs_table)) { # tf_matrix
  tf_matrix[prioritized_tfs_table$tf[i], prioritized_tfs_table$cell_type[i]] <- 1
}

prioritized_snps <- unique(prioritized_table$effect_snp)
snp_matrix <- matrix(0, nrow = length(prioritized_snps), ncol = length(prioritized_cells))
rownames(snp_matrix) <- prioritized_snps
colnames(snp_matrix) <- prioritized_cells

for (i in 1:nrow(prioritized_table)) { # snpp_matrix
  snp_matrix[prioritized_table$effect_snp[i], prioritized_table$cell_type[i]] <- 1
}

prioritized_loci_table <- read.table("prioritized_table.txt", header = TRUE, sep = "\t")
prioritized_loci_table$locus <- paste0("locus_", prioritized_loci_table$locus)
prioritized_loci <- unique(prioritized_loci_table$locus)
loci_matrix <- matrix(0, nrow = length(prioritized_loci), ncol = length(prioritized_cells))
rownames(loci_matrix) <- prioritized_loci
colnames(loci_matrix) <- prioritized_cells

for (i in 1:nrow(prioritized_loci_table)) { # loci_matrix
  loci_matrix[prioritized_loci_table$locus[i], prioritized_loci_table$cell_type[i]] <- 1
}

draw_heatmap <- function(matrix, name) {
  ht <- Heatmap(
    as.matrix(matrix), 
    name = name,
    col = c("0" = "white", "1" = "red"),
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 4),
    row_dend_width = unit(3, "cm"),
    column_dend_height = unit(2, "cm"),
    border = TRUE, # Add borders around the entire heatmap
    rect_gp = gpar(col = "black", lwd = 0.1), # Add borders around the tiles
    row_title = name,
    column_title = "Immune Cell Subtypes",
    row_title_gp = gpar(fontsize = 25), # Increase font size of the y-axis title
    column_title_gp = gpar(fontsize = 25), # Increase font size of the x-axis title
    column_names_gp = gpar(fontsize = 15), # Increase font size of the column names
    heatmap_legend_param = list(
      labels_gp = gpar(fontsize = 15)  # Increase labels font size
    )
  )
  return(ht)
}

gene_heatmap <- draw_heatmap(gene_matrix, "Gene")
tf_heatmap <- draw_heatmap(tf_matrix, "TF")
snp_heatmap <- draw_heatmap(snp_matrix, "Effect-SNP")
loci_heatmap <- draw_heatmap(loci_matrix, "Locus")

png("gene_heatmap.png", width = 2000, height = 2600, res = 300)
draw(gene_heatmap)
dev.off()

png("tf_heatmap.png", width = 2000, height = 2600, res = 300)
draw(tf_heatmap)
dev.off()

png("snp_heatmap.png", width = 2000, height = 2600, res = 300)
draw(snp_heatmap)
dev.off()

png("loci_heatmap.png", width = 2000, height = 2600, res = 300)
draw(loci_heatmap)
dev.off()

# boxplot for variants-to-genes maps
complete_table <- read.csv("complete_table.csv")

complete_table[complete_table == "not found"] <- NA
complete_table[complete_table == "not applicable"] <- NA

complete_table <- complete_table %>% filter(!is.na(directly_mapped_gene) | !is.na(cicero_gene) | !is.na(hic_gene) | !is.na(eqtl_gene))

boxplot_data_com <- data.frame()
for (cell_type in unique(complete_table$cell_type)) {
  cell_type_data <- complete_table[complete_table$cell_type == cell_type,]
  cell_type_data <- cell_type_data %>% group_by(locus) %>% summarise(n = n())
  cell_type_data$cell_type <- cell_type
  boxplot_data_com <- rbind(boxplot_data_com, cell_type_data)
}
colnames(boxplot_data_com) <- c("locus", "n_com", "cell_type")

prioritized_table <- read.table("prioritized_table.txt", header = TRUE, sep = "\t")

boxplot_data_pro <- data.frame()
for (cell_type in unique(prioritized_table$cell_type)) {
  cell_type_data <- prioritized_table[prioritized_table$cell_type == cell_type,]
  cell_type_data <- cell_type_data %>% group_by(locus) %>% summarise(n = n())
  cell_type_data$cell_type <- cell_type
  boxplot_data_pro <- rbind(boxplot_data_pro, cell_type_data)
}
colnames(boxplot_data_pro) <- c("locus", "n_pro", "cell_type")

boxplot_data <- data.frame()
for (cell_type in unique(boxplot_data_com$cell_type)) {
  cell_type_data_com <- boxplot_data_com[boxplot_data_com$cell_type == cell_type,]
  cell_type_data_pro <- boxplot_data_pro[boxplot_data_pro$cell_type == cell_type,]
  cell_type_data <- merge(cell_type_data_com, cell_type_data_pro, by = "locus", all = TRUE)
  cell_type_data[is.na(cell_type_data)] <- 0
  boxplot_data <- rbind(boxplot_data, cell_type_data)
}

boxplot_data <- boxplot_data %>% arrange(cell_type, locus) %>% select(locus, n_com, n_pro, cell_type.x)
colnames(boxplot_data) <- c("locus", "n_com", "n_pro", "cell_type")

df_long <- boxplot_data %>%
  pivot_longer(cols = starts_with("n_"), names_to = "type", values_to = "count")

df_long$locus <- as.factor(df_long$locus)

ggplot(df_long, aes(x = locus, y = count, fill = type)) +
  geom_boxplot() +
  theme_classic() +
  labs(
    title = "Box Plot for Loci",
    x = "Locus",
    y = "Count"
  ) +
  theme(legend.position = "bottom")

ggsave("boxplot.png", width = 20, height = 10, units = "in", dpi = 300)

# histone marks 
complete_table <- read.csv("complete_table.csv")

complete_table <- complete_table %>% separate_rows(rmp_hm_label, sep = ", ")

cell_types <- unique(complete_table$cell_type)
rmps <- unique(complete_table$rmp)
rmp_hm_labels <- unique(complete_table$rmp_hm_label)
rmp_hm_labels <- rmp_hm_labels[!is.na(rmp_hm_labels)]

cell_type_matrix <- matrix(0, nrow = length(rmps), ncol = length(rmp_hm_labels), dimnames = list(rmps, rmp_hm_labels))

cell_type_matrices <- list()

for (cell_type in cell_types) {
  cell_type_table <- complete_table[complete_table$cell_type == cell_type, c("rmp", "rmp_hm_label", "cell_type")]
  rmps_in_cell_type <- unique(cell_type_table$rmp)
  cell_type_matrix <- matrix(0, nrow = length(rmps_in_cell_type), ncol = length(rmp_hm_labels), dimnames = list(rmps_in_cell_type, rmp_hm_labels))
  for (rmp in rmps_in_cell_type) {
    rmp_table <- cell_type_table[cell_type_table$rmp == rmp,]
    for (rmp_hm_label in rmp_hm_labels) {
      if (rmp_hm_label %in% rmp_table$rmp_hm_label) {
        cell_type_matrix[rmp, rmp_hm_label] <- 1
      }
    }
  }
  cell_type_matrices[[cell_type]] <- cell_type_matrix
}

cell_type_matrices <- lapply(cell_type_matrices, function(cell_type_matrix) {
  cell_type_matrix[rowSums(cell_type_matrix) > 0,]
})

cell_type_matrices <- cell_type_matrices[sapply(cell_type_matrices, function(cell_type_matrix) nrow(cell_type_matrix) > 0)]

calculate_enrichment <- function(matrix, cell_type_name) {
  results <- data.frame(
    Cell_Type = character(),
    Label = character(),
    P_Value = numeric(),
    FDR = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (label in colnames(matrix)) {
    presence <- sum(matrix[, label] == 1) # presence (1s) and absence (0s) for the label
    absence <- sum(matrix[, label] == 0)
    total_regions <- nrow(matrix)
    count_table <- matrix(c(presence, total_regions - presence, absence, total_regions - absence), nrow = 2) # 2x2 contingency table
    test_result <- fisher.test(count_table) # Fisher's Exact Test
    p_value <- test_result$p.value
    results <- rbind(results, data.frame(
      Cell_Type = cell_type_name,
      Label = label,
      P_Value = p_value,
      stringsAsFactors = FALSE
    ))
  }
  results$FDR <- p.adjust(results$P_Value, method = "BH") # fdr measurment
  return(results)
}

all_results <- lapply(names(cell_type_matrices), function(cell_type) {
  calculate_enrichment(cell_type_matrices[[cell_type]], cell_type)
})

final_results <- bind_rows(all_results)

write.csv(final_results, "enrichment_hm_labels.csv", row.names = FALSE)

library(qqman) # checking inflation in enrichment analysis
pvals <- final_results$P_Value
qq(pvals, main = "Q-Q Plot of P-values from Enrichment Analysis")
dev.copy(png, "qq_plot.png")
dev.off()

final_results <- final_results %>%
  mutate(NegLogAdjPValue = -log10(FDR))

cell_types <- cell_type_matrices %>% names()

final_matrix <- matrix(0, nrow = length(rmp_hm_labels), ncol = length(cell_types), 
                       dimnames = list(rmp_hm_labels, cell_types))

for (cell_type in cell_types) {
  cell_type_results <- final_results[final_results$Cell_Type == cell_type,]
  for (rmp_hm_label in rmp_hm_labels) {
    if (rmp_hm_label %in% cell_type_results$Label) {
      final_matrix[rmp_hm_label, cell_type] <- cell_type_results[cell_type_results$Label == rmp_hm_label,]$NegLogAdjPValue
    }
  }
}

library(ComplexHeatmap)
library(circlize)

sorted_matrix <- final_matrix[order(rownames(final_matrix)), ]

ht <- Heatmap(
  sorted_matrix, 
  col = circlize::colorRamp2(c(0, median(sorted_matrix), max(sorted_matrix)), c("white", "pink", "red")),
  cluster_columns = TRUE,
  cluster_rows = FALSE, # Disable row clustering
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 15), # Adjust font size for readability
  column_names_gp = gpar(fontsize = 20), # Adjust font size for readability
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(2, "cm"),
  border = TRUE, # Add borders around the entire heatmap
  rect_gp = gpar(col = "black", lwd = 0.5), # Add borders around the tiles
  row_title = "ChromHMM Labels",
  column_title = "Immune Cell Subtypes",
  row_title_gp = gpar(fontsize = 20), # Adjust font size of the y-axis title
  column_title_gp = gpar(fontsize = 35), # Adjust font size of the x-axis title
  heatmap_legend_param = list(
    title = "-Log10(FDR)",
    title_gp = gpar(fontsize = 13),
    labels_gp = gpar(fontsize = 10)  # Adjust labels font size
  )
)

png(filename = "hm_heatmap.png", width = 3600, height = 3600, res = 300)
draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()


