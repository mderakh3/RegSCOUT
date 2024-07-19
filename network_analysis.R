setwd("/Users/m.hossein_drn/Documents/aim_1/network_analysis/")

library(dplyr)
library(tidyr)
library(stringr)

# preparing input for network analysis
prioritized_table <- read.table("prioritized_table.txt", header = TRUE, sep = "\t")

prioritized_table <- prioritized_table %>% separate_rows(tf, sep = ", ")
prioritized_table <- prioritized_table %>% separate_rows(gene_score, sep = ", ")

network_table <- prioritized_table %>% dplyr::select(tf, gene_score, cell_type, rmp_ppa)
network_table$tf_change <- str_extract(network_table$tf, "\\((.*?)\\)")
network_table$tf <- str_remove(network_table$tf, "\\((.*?)\\)")

colnames(network_table) <- c("tf", "gene", "cell_type", "risk_score", "tf_change")
network_table$gene_score <- str_extract(network_table$gene, "\\((.*?)\\)")
network_table$gene <- str_remove(network_table$gene, "\\((.*?)\\)")

cell_types <- unique(network_table$cell_type)

cell_data_list <- list()
for (cell in cell_types) {
  cell_data <- network_table %>% filter(cell_type == cell)
  if (nrow(cell_data) == 0) {
    next
  }
  genes <- data.frame(id = unique(c(cell_data$tf, cell_data$gene)))
  genes$gene_score <- ifelse(genes$id %in% cell_data$gene[cell_data$gene_score == "(3)"], 3, 
                             ifelse(genes$id %in% cell_data$gene[cell_data$gene_score == "(2)"], 2, NA))
  genes$risk_score <- ifelse(genes$id %in% cell_data$gene, cell_data$risk_score, NA)
  genes$type <- ifelse(genes$id %in% cell_data$tf, "TF", "Gene")
  genes$cell_type <- cell
  cell_data_list[[cell]] <- genes
}

cell_data_list <- lapply(cell_data_list, function(x) {
  x <- x %>% separate_rows(id, sep = "::")
  x <- x %>% distinct()
  return(x)
})

for (cell in cell_types) {
  write.table(cell_data_list[[cell]], file = paste0("string_input/", cell, ".txt"), sep = "\t", row.names = FALSE)
}

# processing the output of STRING analysis
string_output_files <- list.files("string_output", pattern = "*.tsv", full.names = TRUE)

string_output_data <- lapply(string_output_files, function(x) {
  read.table(x, sep = "\t", header = FALSE)
})

names(string_output_data) <- gsub(".tsv", "", gsub("string_output/", "", string_output_files))

pathway_files <- list.files("string_output/pathways", pattern = "*.tsv", full.names = TRUE)

pathway_data <- lapply(pathway_files, function(x) {
  read.table(x, sep = "\t", header = FALSE)
})

names(pathway_data) <- gsub("_kegg.tsv", "", gsub("string_output/pathways/", "", pathway_files))

supp_tab_data <- data.frame(matrix(ncol = 5, nrow = 13))
colnames(supp_tab_data) <- c("cell_types", "genes/tfs", "string(score>0.4)", "no_cluster", "major_cluster_stringest_kegg_top5")
supp_tab_data$cell_types <- cell_types

for (cell in cell_types) {
  supp_tab_data[supp_tab_data$cell_types == cell, "genes/tfs"] <- cell_data_list[[cell]] %>% pull(id) %>% unique() %>% length()
  supp_tab_data[supp_tab_data$cell_types == cell, "string(score>0.4)"] <- string_output_data[[cell]] %>% pull(V5) %>% unique() %>% length()
  supp_tab_data[supp_tab_data$cell_types == cell, "no_cluster"] <- string_output_data[[cell]] %>% pull(V2) %>% unique() %>% length()
}

for (cell in cell_types) {
  pathway_data[[cell]]$V5 <- as.numeric(as.character(pathway_data[[cell]]$V5))
  pathway_data[[cell]] <- pathway_data[[cell]] %>% arrange(desc(V5))
  supp_tab_data[supp_tab_data$cell_types == cell, "major_cluster_stringest_kegg"] <- pathway_data[[cell]][c(1:5), "V2"] %>% paste(collapse = ", ")
}

write.table(supp_tab_data, file = "supplementary_table.txt", sep = "\t", row.names = FALSE)
