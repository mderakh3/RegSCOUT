setwd("/Users/m.hossein_drn/Documents/aim_1/prioritizing_maps/")

library(dplyr)
library(stringr)
library(tidyr)

# prioritizing maps 
final_table <- read.csv("complete_table.csv", header = TRUE, sep = ",")
final_table$prioritization_score <- NA
final_table$gene <- NA
  
for (i in 1:nrow(final_table)) {
  final_table$gene[i] <- paste0(final_table$directly_mapped_gene[i], ", ", final_table$cicero_gene[i], ", ", final_table$hic_gene[i], ", ", final_table$eqtl_gene[i])
}

final_table <- final_table %>% separate_rows(gene, sep = ", ")
final_table <- final_table[!final_table$gene %in% "NA",] # remove the rows with "NA" gene

not_found <- final_table[final_table$gene == "not found",]
not_applicable <- final_table[final_table$gene == "not applicable",]
final_table <- final_table[!final_table$gene %in% c("not found", "not applicable"),]
final_table <- distinct(final_table)
rm(not_found, not_applicable)

gene_ppa_df <- final_table %>%
  distinct(cell_type, gene, rmp, .keep_all = TRUE) %>%
  group_by(cell_type, gene) %>%
  summarise(gene_ppa = sum(rmp_ppa, na.rm = TRUE), .groups = 'drop')

final_table <- final_table %>%
  left_join(gene_ppa_df, by = c("cell_type", "gene"))

for (i in 1:nrow(final_table)) {
  gene <- final_table$gene[i]
  cicero_gene <- final_table$cicero_gene[i]
  hic_gene <- final_table$hic_gene[i]
  eqtl_gene <- final_table$eqtl_gene[i]
  if (gene %in% unlist(strsplit(cicero_gene, ", ")) & 
      gene %in% unlist(strsplit(hic_gene, ", ")) & 
      gene %in% unlist(strsplit(eqtl_gene, ", "))) {
    final_table$prioritization_score[i] <- 3
  }
}

nrow(final_table[final_table$prioritization_score %in% 3,])

for (i in 1:nrow(final_table)) {
  gene <- final_table$gene[i]
  directly_mapped_gene <- final_table$directly_mapped_gene[i]
  hic_gene <- final_table$hic_gene[i]
  eqtl_gene <- final_table$eqtl_gene[i]
  if (gene %in% unlist(strsplit(directly_mapped_gene, ", ")) & 
      gene %in% unlist(strsplit(hic_gene, ", ")) & 
      gene %in% unlist(strsplit(eqtl_gene, ", "))) {
    final_table$prioritization_score[i] <- 3
  }
}

nrow(final_table[final_table$prioritization_score %in% 3,])

for (i in 1:nrow(final_table)) { # NK cells don't have gene in intHiC
  gene <- final_table$gene[i]
  cicero_gene <- final_table$cicero_gene[i]
  eqtl_gene <- final_table$eqtl_gene[i]
  hic_interaction_proximity <- final_table$hic_interaction_proximity[i]
  if (final_table$cell_type[i] %in% c("adaptive_NK", "cyto_nk") & 
      gene %in% unlist(strsplit(cicero_gene, ", ")) & 
      gene %in% unlist(strsplit(eqtl_gene, ", ")) & 
      hic_interaction_proximity == "yes") {
    final_table$prioritization_score[i] <- 3
  }
}

nrow(final_table[final_table$prioritization_score %in% 3,])

final_table_3 <- final_table[final_table$prioritization_score %in% 3,]
final_table <- final_table[!final_table$prioritization_score %in% 3,]

for (i in 1:nrow(final_table)) {
  gene <- final_table$gene[i]
  directly_mapped_gene <- final_table$directly_mapped_gene[i]
  hic_gene <- final_table$hic_gene[i]
  if (gene %in% unlist(strsplit(directly_mapped_gene, ", ")) &
      gene %in% unlist(strsplit(hic_gene, ", "))) {
    final_table$prioritization_score[i] <- 2
  }
}

nrow(final_table[final_table$prioritization_score %in% 2,])

for (i in 1:nrow(final_table)) {
  gene <- final_table$gene[i]
  directly_mapped_gene <- final_table$directly_mapped_gene[i]
  eqtl_gene <- final_table$eqtl_gene[i]
  if (gene %in% unlist(strsplit(directly_mapped_gene, ", ")) &
      gene %in% unlist(strsplit(eqtl_gene, ", "))) {
    final_table$prioritization_score[i] <- 2
  }
}

nrow(final_table[final_table$prioritization_score %in% 2,])

final_table_2_1 <- final_table[final_table$prioritization_score %in% 2,]
final_table <- final_table[!final_table$prioritization_score %in% 2,]

for (i in 1:nrow(final_table)) { # NK cells don't have gene in intHiC
  gene <- final_table$gene[i]
  cicero_gene <- final_table$cicero_gene[i]
  hic_interaction_proximity <- final_table$hic_interaction_proximity[i]
  if (final_table$cell_type[i] %in% c("adaptive_NK", "cyto_nk") & 
      gene %in% unlist(strsplit(cicero_gene, ", ")) & 
      hic_interaction_proximity == "yes") {
    final_table$prioritization_score[i] <- 2
  }
}

nrow(final_table[final_table$prioritization_score %in% 2,])

final_table_2_2 <- final_table[final_table$prioritization_score %in% 2,]
final_table <- final_table[!final_table$prioritization_score %in% 2,]

for (i in 1:nrow(final_table)) {
  gene <- final_table$gene[i]
  cicero_gene <- final_table$cicero_gene[i]
  hic_gene <- final_table$hic_gene[i]
  if (gene %in% unlist(strsplit(cicero_gene, ", ")) &
      gene %in% unlist(strsplit(hic_gene, ", "))) {
    final_table$prioritization_score[i] <- 2
  }
}

nrow(final_table[final_table$prioritization_score %in% 2,])

for (i in 1:nrow(final_table)) {
  gene <- final_table$gene[i]
  cicero_gene <- final_table$cicero_gene[i]
  eqtl_gene <- final_table$eqtl_gene[i]
  if (gene %in% unlist(strsplit(cicero_gene, ", ")) &
      gene %in% unlist(strsplit(eqtl_gene, ", "))) {
    final_table$prioritization_score[i] <- 2
  }
}

nrow(final_table[final_table$prioritization_score %in% 2,])

for (i in 1:nrow(final_table)) {
  gene <- final_table$gene[i]
  hic_gene <- final_table$hic_gene[i]
  eqtl_gene <- final_table$eqtl_gene[i]
  if (gene %in% unlist(strsplit(hic_gene, ", ")) &
      gene %in% unlist(strsplit(eqtl_gene, ", "))) {
    final_table$prioritization_score[i] <- 2
  }
}

nrow(final_table[final_table$prioritization_score %in% 2,])

final_table_2_3 <- final_table[final_table$prioritization_score %in% 2,]
final_table <- final_table[!final_table$prioritization_score %in% 2,]

final_table$prioritization_score[is.na(final_table$prioritization_score)] <- 1

prioritized_table <- rbind(final_table, final_table_2_1, final_table_2_2, final_table_2_3, final_table_3)
prioritized_table <- prioritized_table %>% distinct()

write.table(prioritized_table, "scored_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

prioritized_table <- prioritized_table[prioritized_table$gene_ppa > 0.05, , drop = FALSE]
prioritized_table <- prioritized_table[prioritized_table$prioritization_score %in% c(2, 3), , drop = FALSE]

prioritized_table$gene_score <- paste0(prioritized_table$gene, " (", prioritized_table$prioritization_score, ")")

com_score <- prioritized_table %>%
  dplyr::select(-prioritization_score, -cicero_interacting_prom) %>%
  group_by(effect_snp, rmp, cell_type) %>%
  summarise(gene_score = paste(unique(gene_score), collapse = ", "), .groups = "drop")

prioritized_table$gene_score <- NULL
prioritized_table$prioritization_score <- NULL
prioritized_table$gene_ppa <- NULL

prioritized_table <- prioritized_table %>%
  inner_join(com_score, by = c("effect_snp", "rmp", "cell_type")) %>%
  dplyr::select(-gene) %>%
  distinct()

prioritized_table <- prioritized_table %>% separate_rows(directly_mapped_gene, sep = ", ")
prioritized_table <- prioritized_table %>% separate_rows(cicero_gene, sep = ", ")
prioritized_table <- prioritized_table %>% separate_rows(hic_gene, sep = ", ")
prioritized_table <- prioritized_table %>% separate_rows(eqtl_gene, sep = ", ")

processed_data <- prioritized_table %>%
  group_by(cell_type, rmp) %>%
  summarise(
    directly_mapped_gene = paste(unique(directly_mapped_gene), collapse = ", "),
    cicero_gene = paste(unique(cicero_gene), collapse = ", "),
    hic_gene = paste(unique(hic_gene), collapse = ", "),
    eqtl_gene = paste(unique(eqtl_gene), collapse = ", ")
  )

prioritized_table$cicero_interacting_prom <- NULL
prioritized_table$hic_interaction_proximity <- NULL

prioritized_table$cicero_gene <- NA
prioritized_table$hic_gene <- NA
prioritized_table$eqtl_gene <- NA
prioritized_table$directly_mapped_gene <- NA

for (i in 1:nrow(prioritized_table)) {
  current_cell_type <- prioritized_table$cell_type[i]
  current_rmp <- prioritized_table$rmp[i]
  matching_row <- processed_data %>%
    filter(cell_type == current_cell_type & rmp == current_rmp)
  if (nrow(matching_row) > 0) {
    prioritized_table$directly_mapped_gene[i] <- matching_row$directly_mapped_gene
    prioritized_table$cicero_gene[i] <- matching_row$cicero_gene
    prioritized_table$hic_gene[i] <- matching_row$hic_gene
    prioritized_table$eqtl_gene[i] <- matching_row$eqtl_gene
  }
}

prioritized_table <- prioritized_table %>% distinct()
prioritized_table[prioritized_table == "NA"] <- NA
prioritized_table$gene_ppa <- NULL

write.table(prioritized_table, file = "prioritized_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# function to create distribution tables
create_distribution_table <- function(column_name, data) {
  unique_values <- unique(data[[column_name]])
  distribution_table <- data.frame()

  for (value in unique_values) {
    value_table <- data[data[[column_name]] == value,]
    n_prioritized_maps <- nrow(value_table)
    percentage <- n_prioritized_maps / nrow(data) * 100
    distribution_table <- rbind(distribution_table, data.frame(value = value, n_prioritized_maps = n_prioritized_maps, percentage = percentage))
  }

  names(distribution_table)[1] <- column_name
  return(distribution_table)
}

# create distribution tables for cell type, locus, effect snp, and rmp
cell_table_final <- create_distribution_table("cell_type", prioritized_table)
locus_table_final <- create_distribution_table("locus", prioritized_table)
effect_snp_table_final <- create_distribution_table("effect_snp", prioritized_table)
rmp_table_final <- create_distribution_table("rmp", prioritized_table)

write.table(cell_table_final, file = "cell_table_final.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(locus_table_final, file = "locus_table_final.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(effect_snp_table_final, file = "effect_snp_table_final.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rmp_table_final, file = "rmp_table_final.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# create gene distribution table
extract_genes <- function(gene_score) {
  # remove the score in parentheses and split by comma
  genes <- str_remove_all(gene_score, "\\s*\\(.*?\\)") %>% 
    str_split(",\\s*") %>% 
    unlist()
  return(genes)
}

library(purrr)
all_genes <- prioritized_table %>% 
  pull(gene_score) %>% 
  map(extract_genes) %>% 
  unlist()

gene_counts <- tibble(gene = all_genes) %>%
  count(gene, name = "count") %>%
  arrange(desc(count))

write.table(gene_counts, file = "gene_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# create tf distribution table
extract_tfs <- function(tf) {
  # Split by comma and remove any trailing spaces
  tfs <- str_split(tf, ",\\s*") %>% 
    unlist()
  return(tfs)
}

all_tfs <- prioritized_table %>% 
  pull(tf) %>% 
  map(extract_tfs) %>% 
  unlist()

tf_counts <- tibble(tf = all_tfs) %>%
  count(tf, name = "count") %>%
  arrange(desc(count))

tf_counts$tf <- gsub("\\(.*?\\)", "", tf_counts$tf)

tf_counts <- tf_counts %>%
  group_by(tf) %>%
  summarise(count = sum(count)) %>%
  arrange(desc(count))

write.table(tf_counts, file = "tf_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# a table for score = 3 genes in each cell type
score_3_genes <- final_table_3 %>%
  filter(prioritization_score == 3) %>%
  dplyr::select(cell_type, gene) %>% distinct()

n_rows <- score_3_genes %>%
  group_by(gene) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

for (i in 1:nrow(score_3_genes)) {
  current_gene <- score_3_genes$gene[i]
  matching_row <- n_rows %>%
    filter(gene == current_gene)
  score_3_genes$n[i] <- matching_row$n
}

for (i in 1:nrow(score_3_genes)) {
  current_gene <- score_3_genes$gene[i]
  current_n <- score_3_genes$n[i]
  score_3_genes$gene[i] <- paste(current_gene, " (",current_n, ")", sep = "")
}

score_3_genes <- score_3_genes %>%
  group_by(cell_type) %>%
  summarise(gene = paste(unique(gene), collapse = ", "))

write.table(score_3_genes, file = "score_3_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# cell-type specific loci-ppa table
locus_ppa_p_table <- prioritized_table %>%
  dplyr::select(cell_type, locus, rmp, rmp_ppa) %>% distinct()

locus_ppa_sums <- locus_ppa_p_table %>%
  group_by(cell_type, locus) %>%
  summarise(sum_ppa = sum(rmp_ppa, na.rm = TRUE), .groups = 'drop')

max_locus_ppa_sums <- locus_ppa_sums %>%
  group_by(cell_type) %>%
  filter(sum_ppa == max(sum_ppa)) %>%
  ungroup()

write.table(max_locus_ppa_sums, file = "max_locus_ppa_sums.txt", sep = "\t", quote = FALSE, row.names = FALSE)
