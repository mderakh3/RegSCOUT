setwd("/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage/")

library(readr)
library(openxlsx)
library(dplyr)
library(tibble)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)

# preparing input data for eQTL analysis
peak_cluster_matrix_dir <- file.path("/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage", "peak_cluster_matrix.txt")
rmp_snps_dir <- file.path("/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage", "risk_regions_ratio.xlsx")
ci_snps_dir <- file.path("/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage", "ci_snp.txt")
chain_dir = "/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage/hg19ToHg38.over.chain"

eqtl_rmp_cluster_matrix <- read.table(peak_cluster_matrix_dir, sep = "\t", header = TRUE)
eqtl_rmp_cluster_matrix <- rownames_to_column(eqtl_rmp_cluster_matrix, var = "rmp")

rmp_snps <- read.xlsx(rmp_snps_dir, sheet = 1)

rmp_snps <- rmp_snps %>%
  mutate(SNP = gsub(".*-", "", TFSNP)) %>%
  dplyr::select(region, SNP)

ci_snps <- read.delim(ci_snps_dir, header = TRUE)
ci_snps = ci_snps[ci_snps$id %in% rmp_snps$SNP,]

# lifting over hg19 to hg38
chain_data = import.chain(chain_dir) 

ci_snps_granges = GRanges(seqnames = ci_snps$chr, ranges = IRanges(start = ci_snps$pos - 1, end = ci_snps$pos))
ci_snps_granges_hg38 <- liftOver(ci_snps_granges, chain_data)

lifted_snps_df <- data.frame(
  id = ci_snps$id,
  chr = seqnames(ci_snps_granges_hg38),
  start = start(ci_snps_granges_hg38),
  end = end(ci_snps_granges_hg38)
)

for (i in 1:nrow(lifted_snps_df)) {
  lifted_snps_df$allelic_id[i] <- paste(lifted_snps_df$chr.value[i], ":", lifted_snps_df$start.value[i], 
                                        "-", lifted_snps_df$end.value[i], sep = "")
}

# creating allelic_ids
lifted_snps_df$allelic_id <- gsub("chr", "", lifted_snps_df$allelic_id)
lifted_snps_df <- lifted_snps_df %>% 
  dplyr::select(id, allelic_id) %>% 
  merge(rmp_snps, by.x = "id", by.y = "SNP")

eqtl_rmp_cluster_matrix <- eqtl_rmp_cluster_matrix %>%
  merge(lifted_snps_df, by.x = "rmp", by.y = "region") %>%
  dplyr::select(-id)

cell_subtypes <- list(
  cMono = "cMono",
  cd16_Mono = c("ncMono", "iMono"),
  mem_b = "mem_b",
  naive_b = "naive_b",
  act_cd4_t = "act_cd4_t",
  tReg = "tReg",
  cyto_cd8_t = "cyto_cd8_t",
  mem_cd8_t = "mem_cd8_t",
  naive_cd4_t = "naive_cd4_t",
  naive_cd8_t = "naive_cd8_t",
  nk = c("cyto_nk", "adaptive_NK")
)

eqtl_rmp_list <- list()
eqtl_snps_list <- list()

for (subtype in names(cell_subtypes)) {
  columns <- c("rmp", cell_subtypes[[subtype]], "allelic_id")
  eqtl_rmp <- eqtl_rmp_cluster_matrix[, which(colnames(eqtl_rmp_cluster_matrix) %in% columns)]
  
  if (length(cell_subtypes[[subtype]]) > 1) {
    eqtl_rmp <- eqtl_rmp[!(rowSums(eqtl_rmp[, cell_subtypes[[subtype]]]) == 0), ]
  } else {
    eqtl_rmp <- eqtl_rmp[!(eqtl_rmp[, cell_subtypes[[subtype]]] == 0), ]
  }
  
  snps <- as.vector(eqtl_rmp$allelic_id)
  
  eqtl_rmp_list[[subtype]] <- eqtl_rmp
  eqtl_snps_list[[subtype]] <- snps
}

output_dir <- "eQTL_input"
for (subtype in names(eqtl_snps_list)) {
  file_path <- file.path(output_dir, paste0(subtype, "_eqtl_snps.txt"))
  write.table(eqtl_snps_list[[subtype]], file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

print("Now take the eqtl inputs and run them through the bash script to get the eQTL results")

# processing the eQTL results
main_dir <- "/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage/eQTL_output"
header_dir = "/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage/headers.txt"

files <- list.files(path = main_dir, pattern = "*.txt", full.names = TRUE)
file_names <- sub("\\.txt$", "", basename(files))
eqtl_results <- setNames(lapply(files, read.table, header = TRUE), file_names)

header = read.table(header_dir, header = FALSE, sep = "\t")
header$V19 = "allelic_id"
header = header[,c(1:5,19,6:18)]

for (i in 1:length(eqtl_results)) {
  colnames(eqtl_results[[i]]) = header[1,]
}

# combining naive and memory tReg
eqtl_results$tReg_eqtl <- rbind(eqtl_results$tReg_eqtl_nai, eqtl_results$tReg_eqtl_mem) 
eqtl_results$tReg_eqtl_nai <- NULL
eqtl_results$tReg_eqtl_mem <- NULL

# filtering out the eQTL results with p-value > 0.001
for (i in 1:length(eqtl_results)) { 
  eqtl_results[[i]] = eqtl_results[[i]][eqtl_results[[i]]$pvalue < 0.001,]
}

# converting gene_ids to gene symbols
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") 

gene_ids = c()
for (i in 1:length(eqtl_results)) {
  gene_ids = c(gene_ids, eqtl_results[[i]]$gene_id)
}

gene_ids = unique(gene_ids)
attributes_ensemble = listAttributes(ensembl)

gene_names = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = ensembl)

for (i in 1:length(eqtl_results)) {
  eqtl_results[[i]] = merge(eqtl_results[[i]], gene_names, by.x = "gene_id", by.y = "ensembl_gene_id")
}

for (i in 1:length(eqtl_results)) {
  eqtl_results[[i]] = eqtl_results[[i]][,c("gene_id","external_gene_name", "pvalue", "rsid")]
}

output_dir = "eQTL_processed"

for (i in 1:length(eqtl_results)) {
  file_path = file.path(output_dir, paste0(names(eqtl_results)[i], "_processed_results.txt"))
  write.table(eqtl_results[[i]], file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

print("eQTL Done!")
