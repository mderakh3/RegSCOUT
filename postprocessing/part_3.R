setwd("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_3")
library(dplyr)
library(tidyr)
library(stringr)
library(disgenet2r)

all_disgennet_genes <- readRDS("all_disgennet_genes.rds")

disgennet_df <- unique(all_disgennet_genes@qresult[  ,c("gene_symbol", "disease_name","score", "chemicalsIncludedInEvidence", "yearFinal", "diseaseUMLSCUI")] )

disgennet_df$chemicals <- NA
for (i in 1:nrow(disgennet_df)) {
  if (length(disgennet_df$chemicalsIncludedInEvidence[[i]]) > 0) {
    disgennet_df$chemicals[i] <- paste(disgennet_df$chemicalsIncludedInEvidence[[i]], collapse = ", ")
  }
}

enrichr_disease_results <- read.csv("enrichr_disease_results.csv")
enrichr_pathway_results <- read.csv("enrichr_pathway_results.csv")

gene_paths <- enrichr_pathway_results %>%
  separate_rows(Genes, sep = ", ") %>%
  filter(!is.na(Genes)) %>% select(Genes) %>% distinct()

gene_diseases <- enrichr_disease_results %>%
  separate_rows(Genes, sep = ", ") %>%
  filter(!is.na(Genes)) %>% select(Genes) %>% distinct()

all_genes <- gene_paths %>% bind_rows(gene_diseases) %>% distinct()

prioritized_table <- read.table("prioritized_table.txt", header = TRUE, sep = "\t")

prioritized_table <- prioritized_table %>% separate_rows(gene_score, sep = ", ")
prioritized_table$gene_score <- str_remove(prioritized_table$gene_score, "\\(.*\\)")
prioritized_table$gene_score <- str_remove(prioritized_table$gene_score, " ")

prioritized_table <- prioritized_table %>% separate_rows(tf, sep = ", ")
prioritized_table <- prioritized_table %>% separate_rows(tf, sep = "::")
prioritized_table$tf <- str_remove(prioritized_table$tf, "\\(.*\\)")

genes_table = data.frame(Attributes = all_genes$Genes, Type = NA, Cell_type = NA,
                         Pathway = NA, Chemicals = NA, Disease = NA, Loci = NA)

for (i in 1:nrow(genes_table)) {
  if (genes_table$Attributes[i] %in% prioritized_table$tf) {
    genes_table$Type[i] <- "TF"
  }
  if (genes_table$Attributes[i] %in% prioritized_table$gene_score) {
    if (is.na(genes_table$Type[i])) {
      genes_table$Type[i] <- "Gene"
    } else {
      genes_table$Type[i] <- "TF/Gene"
    }
  }
}

for (i in 1:nrow(genes_table)) {
  cell_types <- unique(prioritized_table %>% filter(tf == genes_table$Attributes[i] | gene_score == genes_table$Attributes[i]) %>% select(cell_type) %>% unlist())
  if (length(cell_types) > 0) {
    genes_table$Cell_type[i] <- paste(cell_types, collapse = ", ")
  }
}

enrichr_pathway_results <- enrichr_pathway_results  %>% separate_rows(Genes, sep = ", ")
enrichr_disease_results <- enrichr_disease_results  %>% separate_rows(Genes, sep = ", ")

for (i in 1:nrow(genes_table)) {
  pathways <- unique(enrichr_pathway_results %>% filter(Genes == genes_table$Attributes[i]) %>% select(Pathways) %>% unlist())
  if (length(pathways) > 0) {
    genes_table$Pathway[i] <- paste(pathways, collapse = ", ")
  }
}

for (i in 1:nrow(genes_table)) {
  diseases <- unique(enrichr_disease_results %>% filter(Genes == genes_table$Attributes[i]) %>% select(Disease) %>% unlist())
  if (length(diseases) > 0) {
    genes_table$Disease[i] <- paste(diseases, collapse = ", ")
  }
}

for (i in 1:nrow(genes_table)) {
  chemicals <- unique(disgennet_df %>% filter(gene_symbol == genes_table$Attributes[i]) %>% select(chemicals) %>% unlist())
  if (length(chemicals) > 0) {
    genes_table$Chemicals[i] <- paste(chemicals, collapse = ", ")
  }
}

for (i in 1:nrow(genes_table)) {
  if (genes_table$Type[i] == "Gene" | genes_table$Type[i] == "TF/Gene") {
    loci <- unique(prioritized_table %>% filter(gene_score == genes_table$Attributes[i]) %>% select(locus) %>% unlist())
    if (length(loci) > 0) {
      genes_table$Loci[i] <- paste(loci, collapse = ", ")
    }
  }
}


genes_table[genes_table == "NA"] <- NA

library(openxlsx)
write.xlsx(genes_table, "genes_table.xlsx", rowNames = FALSE)

# loci based table 
loci_based_table <- data.frame(matrix(ncol = 4, nrow = 33))
colnames(loci_based_table) <- c("Loci", "Genes", "Pathways", "Diseases")

loci_based_table$Loci <- unique(prioritized_table$locus)

processed_data <- prioritized_table %>%
  group_by(locus) %>%
  summarize(
    gene_score = paste(unique(gene_score), collapse = ", "),
    tf = paste(unique(tf), collapse = ", ")
  )

for (i in 1:nrow(processed_data)) {
  if (processed_data$gene_score[i] == "") {
    processed_data$gene_score[i] <- processed_data$tf[i]
  } else {
    processed_data$gene_score[i] <- paste(processed_data$gene_score[i], processed_data$tf[i], sep = ", ")
  }
}

for (i in 1:nrow(loci_based_table)) {
  genes <- unique(processed_data %>% filter(locus == loci_based_table$Loci[i]) %>% select(gene_score) %>% unlist())
  if (length(genes) > 0) {
    loci_based_table$Genes[i] <- paste(genes, collapse = ", ")
  }
}

for (i in 1:nrow(loci_based_table)) {
  pathways <- unique(enrichr_pathway_results %>% filter(Genes %in% unlist(strsplit(loci_based_table$Genes[i], ", "))) %>% select(Pathways) %>% unlist())
  if (length(pathways) > 0) {
    loci_based_table$Pathways[i] <- paste(pathways, collapse = ", ")
  }
}

for (i in 1:nrow(loci_based_table)) {
  diseases <- unique(enrichr_disease_results %>% filter(Genes %in% unlist(strsplit(loci_based_table$Genes[i], ", "))) %>% select(Disease) %>% unlist())
  if (length(diseases) > 0) {
    loci_based_table$Diseases[i] <- paste(diseases, collapse = ", ")
  }
}

loci_based_table <- loci_based_table[order(-table(loci_based_table$Loci)),]

write.xlsx(loci_based_table, "loci_based_table.xlsx", rowNames = FALSE)

# pathway based table
pathway_based_table <- data.frame(matrix(ncol = 3, nrow = 10))
colnames(pathway_based_table) <- c("Pathways", "Genes", "Loci")

pathway_based_table$Pathways <- unique(enrichr_pathway_results$Pathways)

for (i in 1:nrow(pathway_based_table)) {
  genes <- unique(enrichr_pathway_results %>% filter(Pathways == pathway_based_table$Pathways[i]) %>% select(Genes) %>% unlist())
  if (length(genes) > 0) {
    pathway_based_table$Genes[i] <- paste(genes, collapse = ", ")
  }
}

for (i in 1:nrow(pathway_based_table)) {
  loci <- unique(prioritized_table %>% filter(gene_score %in% unlist(strsplit(pathway_based_table$Genes[i], ", "))) %>% select(locus) %>% unlist())
  if (length(loci) > 0) {
    pathway_based_table$Loci[i] <- paste(loci, collapse = ", ")
  }
}

for (i in 1:nrow(pathway_based_table)) {
  if (is.na(pathway_based_table$Loci[i])) {
    loci <- unique(prioritized_table %>% filter(tf %in% unlist(strsplit(pathway_based_table$Genes[i], ", "))) %>% select(locus) %>% unlist())
    pathway_based_table$Loci[i] <- paste(loci, collapse = ", ")
  }
}

write.xlsx(pathway_based_table, "pathway_based_table.xlsx", rowNames = FALSE)
