setwd("/Users/m.hossein_drn/Documents/aim_1/gene_set_enrichment/")

library(readxl)
library(dplyr)
library(tidyverse)
library(xlsx)

# preparing input for gene set enrichment analysis
prioritized_table <- read.table("prioritized_table.txt", header = TRUE, sep = "\t")

prioritized_table <- prioritized_table %>% separate_rows(tf, sep = ", ")
prioritized_table <- prioritized_table %>% separate_rows(gene_score, sep = ", ")

genes_table <- prioritized_table %>% dplyr::select(tf, gene_score, cell_type)
genes_table$tf <- str_remove(genes_table$tf, "\\((.*?)\\)")
genes_table <- genes_table %>% separate_rows(tf, sep = "::")
genes_table$gene_score <- str_remove(genes_table$gene_score, "\\((.*?)\\)")

cell_types <- unique(genes_table$cell_type)

gene.lists <- list()
for (cell in cell_types) {
  cell_data <- genes_table %>% filter(cell_type == cell)
  if (nrow(cell_data) == 0) {
    next
  }
  genes <- unique(c(cell_data$tf, cell_data$gene_score))
  gene.lists[[cell]] <- genes
}

gene.lists <- lapply(gene.lists, function(x) {
  x <- x %>% str_replace_all(" ", "")
  return(x)
})

# function for gene set enrichment analysis
require(httr)
require(rjson)
require(dplyr)
require(imager)

register.gene.list <- function(gene.list)
{
  print("Registering the gene list...")
  Sys.sleep(3)
  ENRICHR_URL <- 'https://maayanlab.cloud/Enrichr/addList'
  gene.list.str <- Reduce(function(x, y){paste(x, y, sep='\n')}, gene.list)
  response <- POST(ENRICHR_URL, body = list(list = gene.list.str, description = ''))
  #print(response)
  fromJSON(content(response, 'text'))$userListId
}

enrichr <- function(user.list.id, gene_set_library=c('MSigDB_Hallmark_2020', 'DisGeNET'))
{
  ENRICHR_URL <- 'https://maayanlab.cloud/Enrichr/enrich'
  response <- GET(
    ENRICHR_URL,
    query=list(userListId = user.list.id, backgroundType = gene_set_library)
  )
  enriched.res <- fromJSON(content(response, 'text'))
  enriched.res <- Reduce(rbind, enriched.res[[1]]) %>% as.data.frame()
  if(length(enriched.res) == 0)
  {
    enriched.res <- data.frame(Name=double(), P.value=double(), Odds.ratio=double(), Combined.score=double(), Genes=double(), Adj.p.val=double())
  }else{
    enriched.res <- enriched.res[,2:7]
    colnames(enriched.res) <- c('Name', 'P.value', 'Odds.ratio', 'Combined.score', 'Genes', 'Adj.p.val')
    rownames(enriched.res) <- NULL
  }
  enriched.res
}

# pathways and diseases
enrichr.path.result <- gene.lists %>% lapply(function(gene.list){
  gene.list.id <- register.gene.list(gene.list)
  enrichr.res <- enrichr(gene.list.id, gene_set_library= 'MSigDB_Hallmark_2020')
  enrichr.res <- enrichr.res %>% dplyr::select(Name, Genes, Adj.p.val) %>% filter(Adj.p.val < .05)
  enrichr.res
})

saveRDS(enrichr.path.result, "enrichr_pathway_results.rds")

enrichr.dis.result <- gene.lists %>% lapply(function(gene.list){
  gene.list.id <- register.gene.list(gene.list)
  enrichr.res <- enrichr(gene.list.id, gene_set_library= 'DisGeNET')
  enrichr.res <- enrichr.res %>% dplyr::select(Name, Genes, Adj.p.val) %>% filter(Adj.p.val < .05)
  enrichr.res
})

saveRDS(enrichr.dis.result, "enrichr_disease_results.rds")

dis_filt_result <- lapply(enrichr.dis.result, function(df) {
  df <- df %>% filter(grepl("Inflammatory Bowel Diseases|Crohn|Colitis", Name))
  df
})

results <- data.frame(Cell_Type = character(),
                      Disease = character(),
                      Adjusted_P_Value = numeric(),
                      Genes = character(),
                      stringsAsFactors = FALSE)

for (cell_type in names(dis_filt_result)) {
  
  cell_data <- dis_filt_result[[cell_type]]
  
  if (length(cell_data$Name) == 0) next
  
  diseases <- cell_data$Name
  adjusted_p_values <- cell_data$Adj.p.val
  genes <- sapply(cell_data$Genes, paste, collapse = ", ")
  
  cell_type_df <- data.frame(Cell_Type = cell_type,
                             Disease = unlist(diseases),
                             Adjusted_P_Value = unlist(adjusted_p_values),
                             Genes = genes,
                             stringsAsFactors = FALSE)
  
  results <- rbind(results, cell_type_df)
}

results$Gene_Ratio <- sapply(results$Genes, function(x) {
  length(unlist(strsplit(x, ", "))) / nrow(results)
})

write.csv(results, "enrichr_disease_results.csv", row.names = FALSE)

results_paths <- data.frame(Cell_Type = character(),
                      Pathways = character(),
                      Adjusted_P_Value = numeric(),
                      Genes = character(),
                      stringsAsFactors = FALSE)

for (cell_type in names(enrichr.path.result)) {
  
  cell_data <- enrichr.path.result[[cell_type]]
  
  if (length(cell_data$Name) == 0) next
  
  pathways <- cell_data$Name
  adjusted_p_values <- cell_data$Adj.p.val
  genes <- sapply(cell_data$Genes, paste, collapse = ", ")
  
  cell_type_df <- data.frame(Cell_Type = cell_type,
                             Pathways = unlist(pathways),
                             Adjusted_P_Value = unlist(adjusted_p_values),
                             Genes = genes,
                             stringsAsFactors = FALSE)
  
  results_paths <- rbind(results_paths, cell_type_df)
}

results_paths$Gene_Ratio <- sapply(results_paths$Genes, function(x) {
  length(unlist(strsplit(x, ", "))) / nrow(results_paths)
})

write.csv(results_paths, "enrichr_pathway_results.csv", row.names = FALSE)

# bubble plot
results <- results %>%
  mutate(NegLogAdjPValue = -log10(Adjusted_P_Value))

ggplot(results, aes(x = Disease, y = Cell_Type, size = Gene_Ratio, color = NegLogAdjPValue)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Gene Ratio", range = c(7, 20)) +
  scale_color_gradient(low = "pink", high = "red", name = "-log10(Adjusted P-Value)") +
  theme_minimal() +
  labs(
    x = "Diseases",
    y = "Cell Subtypes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )

ggsave("disease_bubble_plot.png", width = 12, height = 12, dpi = 300)

results_paths <- results_paths %>%
  mutate(NegLogAdjPValue = -log10(Adjusted_P_Value))

ggplot(results_paths, aes(x = Pathways, y = Cell_Type, size = Gene_Ratio, color = NegLogAdjPValue)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Gene Ratio", range = c(7, 20)) +
  scale_color_gradient(low = "pink", high = "red", name = "-log10(Adjusted P-Value)") +
  theme_minimal() +
  labs(
    x = "Pathways",
    y = "Cell Subtypes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )

ggsave("pathway_bubble_plot.png", width = 12, height = 12, dpi = 300)

