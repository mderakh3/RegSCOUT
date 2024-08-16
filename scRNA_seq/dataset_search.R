setwd("/Users/m.hossein_drn/Documents/aim_1/scRNA_seq/")

library(rentrez)
library(dplyr)
library(stringr)
library(rvest)
library(dplyr)
library(openxlsx)

search_terms <- "(IBD[TIAB] OR Ulcerative Colitis[TIAB] OR Crohn[TIAB] OR Inflammatory Bowel Disease[TIAB]) AND (single-cell RNA sequencing[TIAB] OR scRNA[TIAB] OR single-cell RNA seq[TIAB] OR single-cell RNA-seq[TIAB])"

search_results <- entrez_search(db = "pubmed", term = search_terms, retmax = 500)

# fetch details for each article
article_details <- lapply(search_results$ids, function(id) {
  article <- entrez_summary(db = "pubmed", id = id)
  title <- article$title
  journal <- article$fulljournalname
  year <- article$pubdate
  keywords <- c("IBD", "Ulcerative Colitis", "Crohn", "Inflammatory Bowel Disease", "single-cell RNA sequencing", "scRNA", "single-cell RNA seq", "single-cell RNA-seq") %>%
    str_c(collapse = "|") %>%
    str_extract_all(article$source)
  
  data.frame(
    title = title,
    journal = journal,
    year = year,
    keywords = paste(unique(unlist(keywords)), collapse = ", ")
  )
})

results_df <- bind_rows(article_details)

write.xlsx(results_df, "scRNA_seq_datasets.xlsx", row.names = FALSE)
