setwd("/Users/m.hossein_drn/Documents/aim_1/final_outputs/")
library(tibble)
library(xlsx)
library(stringr)
library(tidyr)
library(dplyr)
library(readxl)
library(GenomicRanges)

# creating the final output data frame
final_output <- data.frame(matrix(ncol = 18, nrow = 1))

colnames(final_output) <- c("locus",	"trait", "loci_region",	"lead",	"lead_ppa",	
                            "effect_snp",	"effect_ppa",	"tf","rmp",	"rmp_hm_label", "rmp_ppa", 
                            "cicero_interacting_prom", "hic_interaction_proximity",
                            "directly_mapped_gene",	"cicero_gene",	"hic_gene",	
                            "eqtl_gene", "cell_type")

# assigning the locus number to each ci_snp
fine_loci_head <- read.xlsx("/Users/m.hossein_drn/Documents/aim_1/final_outputs/finemapping_outputs/fine_loci_head.xlsx", sheetIndex = 1)
ci_snp <- read.delim("/Users/m.hossein_drn/Documents/aim_1/final_outputs/finemapping_outputs/ci_snp.txt")

fine_loci_head$credibleset <- NULL

loci_regions_grange <- GRanges(seqnames = paste0("chr", fine_loci_head$CHR), 
                                  ranges = IRanges(start = fine_loci_head$credible_start, 
                                                               end = fine_loci_head$credible_end))
loci_regions_grange$locus_id <- fine_loci_head$locus

snp_grange <- GRanges(seqnames = ci_snp$chr, 
                      ranges = IRanges(start = ci_snp$pos - 1, 
                                       end = ci_snp$pos))

overlaps <- findOverlaps(snp_grange, loci_regions_grange)

snp_with_locus <- data.frame(
  snp_index = queryHits(overlaps),
  locus_id = loci_regions_grange$locus_id[subjectHits(overlaps)]
)

snp_grange$locus_id <- NA
snp_grange$locus_id[snp_with_locus$snp_index] <- snp_with_locus$locus_id

snp_grange <- data.frame(snp_grange)

ci_snp$locus <- NA

ci_snp$pos <- as.numeric(ci_snp$pos)
ci_snp$chr <- as.character(ci_snp$chr)

for (i in 1:nrow(ci_snp)) {
  chr <- as.character(ci_snp$chr[i])
  pos <- as.numeric(ci_snp$pos[i])
  match_index <- which(snp_grange$seqnames == chr & snp_grange$end == pos)
  if (length(match_index) > 0) {
    ci_snp$locus[i] <- snp_grange$locus_id[match_index]
  }
}

rm(snp_grange, snp_with_locus, overlaps, loci_regions_grange)
gc()

write.table(ci_snp, file = "/Users/m.hossein_drn/Documents/aim_1/final_outputs/ci_snp_locus.txt", sep = "\t", row.names = FALSE)

# preparing tf data for effect_snps located on rmps
risk_regions_ratio <- read.xlsx("/Users/m.hossein_drn/Documents/aim_1/final_outputs/scATAC_outputs/risk_regions_ratio.xlsx", sheetIndex = 1)
risk_regions_ratio$log_lik_ratio <- NA

log_lik_ratio <- read.delim("/Users/m.hossein_drn/Documents/aim_1/final_outputs/TF_outputs/log_lik_ratio.txt")
log_lik_ratio <- tibble::rownames_to_column(log_lik_ratio, "TFSNP")

# adding the log_lik_ratio column to the risk_regions_ratio data frame
for (i in 1:nrow(risk_regions_ratio)) {
  tfsnp <- risk_regions_ratio$TFSNP[i]
  match_index <- which(log_lik_ratio$TFSNP == tfsnp)
  if (length(match_index) > 0) {
    risk_regions_ratio$log_lik_ratio[i] <- log_lik_ratio$log_lik_ratio[match_index]
  }
}

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

for (i in 1:nrow(risk_regions_ratio)) {
  if (risk_regions_ratio$log_lik_ratio[i] > 0) {
    risk_regions_ratio$direction[i] <- "(+)"
  } else {
    risk_regions_ratio$direction[i] <- "(-)"
  }
}

risk_regions_ratio$tf_direction <- paste0(risk_regions_ratio$tf, risk_regions_ratio$direction)

# addding snps and tfs information to the final_output data frame
final_output <- final_output[rep(1, nrow(risk_regions_ratio)),]
final_output$effect_snp <- risk_regions_ratio$snp

for (i in 1:nrow(risk_regions_ratio)) {
  current_tf <- risk_regions_ratio$tf_direction[i]
  current_snp <- risk_regions_ratio$snp[i]
  match_index <- which(risk_regions_ratio$snp == current_snp & risk_regions_ratio$tf_direction == current_tf)
  if (length(match_index) > 0) {
    final_output$tf[i] <- risk_regions_ratio$tf_direction[match_index]
  }
}

combined_tfs <- final_output %>% group_by(effect_snp) %>% summarise(tf = paste(tf, collapse = ", "))
final_output$tf <- NA

for (i in 1:nrow(final_output)) {
  current_snp <- final_output$effect_snp[i]
  match_index <- which(combined_tfs$effect_snp == current_snp)
  if (length(match_index) > 0) {
    final_output$tf[i] <- combined_tfs$tf[match_index]
  }
}

final_output <- final_output[!duplicated(final_output$effect_snp),]

for (i in 1:nrow(final_output)) {
  current_snp <- final_output$effect_snp[i]
  match_index <- which(ci_snp$id == current_snp)
  if (length(match_index) > 0) {
    final_output$effect_ppa[i] <- ci_snp$PPA[match_index]
    final_output$locus[i] <- ci_snp$locus[match_index]
  }
}

for (i in 1:nrow(fine_loci_head)) {
  fine_loci_head$loci_region[i] <- paste0("chr", fine_loci_head$CHR[i], "-", fine_loci_head$credible_start[i], "-", fine_loci_head$credible_end[i])
}

for (i in 1:nrow(final_output)) {
  current_locus <- final_output$locus[i]
  match_index <- which(fine_loci_head$locus == current_locus)
  if (length(match_index) > 0) {
    final_output$loci_region[i] <- fine_loci_head$loci_region[match_index]
    final_output$lead[i] <- fine_loci_head$LeadSNP[match_index]
  }
}

for (i in 1:nrow(final_output)) {
  current_lead <- final_output$lead[i]
  match_index <- which(ci_snp$id == current_lead)
  if (length(match_index) > 0) {
    final_output$lead_ppa[i] <- ci_snp$PPA[match_index]
  }
}

for (i in 1:nrow(final_output)) {
  current_lead <- final_output$lead[i]
  match_index <- which(ci_snp$id == current_lead)
  if (length(match_index) > 0) {
    final_output$trait[i] <- ci_snp$trait[match_index]
  }
}

final_output <- final_output[order(final_output$locus),]
rownames(final_output) <- NULL

# adding rmp information
snp_rmp <- risk_regions_ratio %>% group_by(snp) %>% summarise(rmp = paste0(region, collapse = ", "))
snp_rmp$rmp <- gsub(",.*", "", snp_rmp$rmp)

for (i in 1:nrow(final_output)) {
  current_snp <- final_output$effect_snp[i]
  match_index <- which(snp_rmp$snp == current_snp)
  if (length(match_index) > 0) {
    final_output$rmp[i] <- snp_rmp$rmp[match_index]
  }
}

risk_regions_ppa <- read.xlsx("/Users/m.hossein_drn/Documents/aim_1/final_outputs/scATAC_outputs/risk_regions_ppa.xlsx", sheetIndex = 1)
risk_regions_ppa <- risk_regions_ppa %>% dplyr::select(region, sumPPA, cell) %>% distinct()

for (i in 1:nrow(final_output)) {
  current_rmp <- final_output$rmp[i]
  match_index <- which(risk_regions_ppa$region == current_rmp)
  if (length(match_index) > 0) {
    final_output$rmp_ppa[i] <- risk_regions_ppa$sumPPA[match_index]
    final_output$cell_type[i] <- risk_regions_ppa$cell[match_index]
  }
}

# adding promoter-linked rmps genes information
rmp_promoter <- read.xlsx("/Users/m.hossein_drn/Documents/aim_1/final_outputs/scATAC_outputs/rmp_promoter_data.xlsx", sheetIndex = 1)

rmp_promoter$prom_gene = paste(rmp_promoter$rmp_prom, rmp_promoter$gene_name, sep = "_")
rmp_promoter <- rmp_promoter[!duplicated(rmp_promoter$prom_gene),]
rmp_promoter <- rmp_promoter[,c("rmp_prom", "gene_name")]

rmp_gene <- rmp_promoter %>% group_by(rmp_prom) %>% summarise(gene = paste0(gene_name, collapse = ", "))

for (i in 1:nrow(final_output)) {
  current_rmp <- final_output$rmp[i]
  match_index <- which(rmp_gene$rmp_prom == current_rmp)
  if (length(match_index) > 0) {
    final_output$directly_mapped_gene[i] <- rmp_gene$gene[match_index]
  } else {
    final_output$directly_mapped_gene[i] <- "-"
  }
}

# adding cicero information
final_output <- final_output %>% separate_rows(cell_type, sep = ",")

cell_subtypes <- final_output$cell_type %>% unique()

cicero <- read.xlsx("/Users/m.hossein_drn/Documents/aim_1/final_outputs/scATAC_outputs/cic_rmp_enhancers.xlsx", sheetIndex = 1)
cicero$genes <- gsub(",", ", ", cicero$genes)

peaks <- "/Users/m.hossein_drn/Documents/aim_1/final_outputs/scATAC_outputs/journal.pgen.1010759.s019.xlsx"
peaks_data <- read_excel(peaks)

peaks_data$regions = paste0(peaks_data$chr,"-", peaks_data$start, "-", peaks_data$end)
peaks_data <- peaks_data %>% separate_rows(`cell sub-types`, sep = ",") # split the 'cell sub-types'

cell_types <- c("t", "b", "nk", "mono")
cell_type_names <- c("t_cic", "b_cic", "nk_cic", "mono_cic")

for (i in seq_along(cell_types)) {
  cell <- cell_types[i]
  cell_name <- cell_type_names[i]
  temp_df <- cicero[, c("Peak1", "Peak2", cell, "genes", "promoter", "enhancer")]
  temp_df <- temp_df[!is.na(temp_df[[cell]]), ]
  assign(cell_name, temp_df)
}

peak_names <- paste0(cell_subtypes, "_peaks")

for (i in seq_along(cell_subtypes)) { # finding peaks in each cell subtype
  cell_type <- cell_subtypes[i]
  peak_name <- peak_names[i]
  temp_peaks <- peaks_data[peaks_data$`cell sub-types` == cell_type, "regions"]
  assign(peak_name, temp_peaks)
}

act_cd4_t_cic = t_cic[t_cic$Peak1 %in% act_cd4_t_peaks$regions & t_cic$Peak2 %in% act_cd4_t_peaks$regions,] 
naive_cd4_t_cic = t_cic[t_cic$Peak1 %in% naive_cd4_t_peaks$regions & t_cic$Peak2 %in% naive_cd4_t_peaks$regions,]
cyto_cd8_t_cic = t_cic[t_cic$Peak1 %in% cyto_cd8_t_peaks$regions & t_cic$Peak2 %in% cyto_cd8_t_peaks$regions,]
naive_cd8_t_cic = t_cic[t_cic$Peak1 %in% naive_cd8_t_peaks$regions & t_cic$Peak2 %in% naive_cd8_t_peaks$regions,]
mem_cd8_t_cic = t_cic[t_cic$Peak1 %in% mem_cd8_t_peaks$regions & t_cic$Peak2 %in% mem_cd8_t_peaks$regions,]
tReg_cic = t_cic[t_cic$Peak1 %in% tReg_peaks$regions & t_cic$Peak2 %in% tReg_peaks$regions,]
ncMono_cic = mono_cic[mono_cic$Peak1 %in% ncMono_peaks$regions & mono_cic$Peak2 %in% ncMono_peaks$regions,]
cMono_cic = mono_cic[mono_cic$Peak1 %in% cMono_peaks$regions & mono_cic$Peak2 %in% cMono_peaks$regions,]
iMono_cic = mono_cic[mono_cic$Peak1 %in% iMono_peaks$regions & mono_cic$Peak2 %in% iMono_peaks$regions,]
adaptive_NK_cic = nk_cic[nk_cic$Peak1 %in% adaptive_NK_peaks$regions & nk_cic$Peak2 %in% adaptive_NK_peaks$regions,]
cyto_nk_cic = nk_cic[nk_cic$Peak1 %in% cyto_nk_peaks$regions & nk_cic$Peak2 %in% cyto_nk_peaks$regions,]
mem_b_cic = b_cic[b_cic$Peak1 %in% mem_b_peaks$regions & b_cic$Peak2 %in% mem_b_peaks$regions,]
naive_b_cic = b_cic[b_cic$Peak1 %in% naive_b_peaks$regions & b_cic$Peak2 %in% naive_b_peaks$regions,]

cic_data_names <- c("act_cd4_t_cic", "naive_cd4_t_cic", "ncMono_cic", "cMono_cic", 
                          "iMono_cic", "adaptive_NK_cic", "cyto_nk_cic", "mem_b_cic", 
                          "naive_b_cic", "cyto_cd8_t_cic", "naive_cd8_t_cic", "mem_cd8_t_cic", 
                          "tReg_cic")

for (i in seq_along(cic_data_names)) {
  original_df_name <- cic_data_names[i]
  original_df <- get(original_df_name)
  selected_df <- original_df[, c("enhancer", "promoter", "genes")]
  assign(original_df_name, selected_df)
}

for (df_name in cic_data_names) {
  df <- get(df_name)
  df <- df %>%
    separate_rows(genes, sep = ", ")
  assign(df_name, df)
}


for (df_name in cic_data_names) {
  df <- get(df_name)
  df$prom_gene <- paste0(df$promoter, "_", df$genes)
  assign(df_name, df)
}

for (df_name in cic_data_names) {
  df <- get(df_name)
  df <- df %>%
    group_by(enhancer) %>%
    summarise(
      prom_gene = paste(unique(prom_gene), collapse = ", "),
      promoter = unique(promoter)[1],
      genes = paste(unique(genes), collapse = ", ")
    )
  df <- as.data.frame(df)
  assign(df_name, df)
}

cells_without_cicero <- c("plasma", "mkc", "cDC", "pDC")

final_output$prom_gene <- NA 

for (i in 1:nrow(final_output)) {
  cell_subtype <- final_output$cell_type[i]
  if (cell_subtype %in% cells_without_cicero) {
    next
  }
  match_cic_name <- paste0(cell_subtype, "_cic")
  match_cic <- get(match_cic_name)
  if (final_output$rmp[i] %in% match_cic$enhancer) {
    final_output$prom_gene[i] <- match_cic$prom_gene[match_cic$enhancer == final_output$rmp[i]]
  }
}

final_output <- final_output %>% separate_rows(prom_gene, sep = ", ") 
final_output <- final_output %>% separate(prom_gene, c("promoter", "cicero_gene"), sep = "_", remove = F)
final_output$prom_gene <- NULL

promoter_genes <- final_output %>% group_by(promoter) %>% summarise(cicero_gene = paste(unique(cicero_gene), collapse = ", "), .groups = "drop")
final_output$cicero_gene <- NA

for (i in 1:nrow(promoter_genes)) {
  promoter <- promoter_genes$promoter[i]
  genes <- promoter_genes$cicero_gene[i]
  final_output$cicero_gene[final_output$promoter == promoter] <- genes
  final_output$cicero_interacting_prom[final_output$promoter == promoter] <- promoter
}

final_output <- final_output %>% distinct()
final_output$final_output <- NULL
final_output$promoter <- NULL
final_output$cicero_gene[is.na(final_output$cicero_gene)] <- "-"

# adding hic data
base_dir <- "/Users/m.hossein_drn/Documents/aim_1/final_outputs/HiC_outputs/"
file_names <- c("b_pcHiC_Cic.txt", "t_pcHiC_Cic.txt", "nk_intHiC_Cic.txt", "mono_pcHiC_Cic.txt")
df_names <- c("b_hic", "t_hic", "nk_hic", "mono_hic")

for (i in seq_along(file_names)) {
  file_path <- paste0(base_dir, file_names[i])
  df_name <- df_names[i]
  df <- read.table(file_path, header = TRUE, sep = "\t")
  if (df_name != "nk_hic") {
    df$genes <- paste(df$pchic_ba_genes, df$pchic_oe_genes, sep = ", ")
    df$genes <- gsub("\\.", "", df$genes)
    df$genes <- gsub(", $", "", df$genes)
    df$genes <- gsub(";", ", ", df$genes)
  }
  assign(df_name, df)
}

mapping <- list(
  act_cd4_t_hic = list(data = "t_hic", peaks = "act_cd4_t_peaks"),
  naive_cd4_t_hic = list(data = "t_hic", peaks = "naive_cd4_t_peaks"),
  ncMono_hic = list(data = "mono_hic", peaks = "ncMono_peaks"),
  cMono_hic = list(data = "mono_hic", peaks = "cMono_peaks"),
  iMono_hic = list(data = "mono_hic", peaks = "iMono_peaks"),
  adaptive_NK_hic = list(data = "nk_hic", peaks = "adaptive_NK_peaks"),
  cyto_nk_hic = list(data = "nk_hic", peaks = "cyto_nk_peaks"),
  mem_b_hic = list(data = "b_hic", peaks = "mem_b_peaks"),
  naive_b_hic = list(data = "b_hic", peaks = "naive_b_peaks"),
  cyto_cd8_t_hic = list(data = "t_hic", peaks = "cyto_cd8_t_peaks"),
  naive_cd8_t_hic = list(data = "t_hic", peaks = "naive_cd8_t_peaks"),
  mem_cd8_t_hic = list(data = "t_hic", peaks = "mem_cd8_t_peaks"),
  tReg_hic = list(data = "t_hic", peaks = "tReg_peaks")
)

for (new_df_name in names(mapping)) {
  data_name <- mapping[[new_df_name]]$data
  peaks_name <- mapping[[new_df_name]]$peaks
  df <- get(data_name)
  peaks <- get(peaks_name)
  filtered_df <- df[df$cic_promoter %in% peaks$regions & df$cic_enhancer %in% peaks$regions, ]
  assign(new_df_name, filtered_df)
}

pchic_df_names <- c("act_cd4_t_hic", "naive_cd4_t_hic", "ncMono_hic", "cMono_hic", 
                  "iMono_hic", "mem_b_hic", 
                  "naive_b_hic", "cyto_cd8_t_hic", "naive_cd8_t_hic", "mem_cd8_t_hic", 
                  "tReg_hic")

for (df_name in pchic_df_names) {
  df <- get(df_name)
  df$loop <- paste0(df$cic_enhancer, "_", df$cic_promoter)
  df <- df[!duplicated(df$loop), ]
  df <- df[, c("cic_enhancer", "cic_promoter", "genes")]
  assign(df_name, df)
}

inthic_df_names <- c("adaptive_NK_hic", "cyto_nk_hic") # because they didn't have genes inintHiC

for (df_name in inthic_df_names) {
  df <- get(df_name)
  df$loop <- paste0(df$cic_enhancer, "_", df$cic_promoter)
  df <- df[!duplicated(df$loop), ]
  df <- df[, c("cic_enhancer", "cic_promoter", "cic_genes")]
  assign(df_name, df)
}

for (df_name in inthic_df_names) {
  df <- get(df_name)
  df$genes <- "not applicable"
  assign(df_name, df)
}

hic_df_names <- c(pchic_df_names, inthic_df_names)

for (df_name in hic_df_names) {
  df <- get(df_name)
  df$loop <- paste0(df$cic_enhancer, "_", df$cic_promoter)
  assign(df_name, df)
}

for (i in 1:nrow(final_output)) {
  if (!is.na(final_output$cicero_interacting_prom[i])) {
    cell_subtype <- final_output$cell_type[i]
    if (cell_subtype %in% cells_without_cicero) {
      next
    }
    match_hic_name <- paste0(cell_subtype, "_hic")
    match_hic <- get(match_hic_name)
    loop <- paste0(final_output$rmp[i], "_", final_output$cicero_interacting_prom[i])
    if (loop %in% match_hic$loop) {
      final_output$hic_gene[i] <- match_hic$genes[match_hic$loop == loop]
      final_output$hic_interaction_proximity[i] <- "yes"
    } else {
      final_output$hic_interaction_proximity[i] <- "no"
    }
  }
}

for (i in 1:nrow(final_output)) {
  if (!is.na(final_output$hic_interaction_proximity[i])) {
    if (final_output$hic_interaction_proximity[i] == "no") {
      if (is.na(final_output$hic_gene[i])) {
        final_output$hic_gene[i] <- "not found"
      }
    }
  }
}

final_output$hic_gene[is.na(final_output$hic_gene)] <- "-"

# adding histone marks data for rmps
base_dir <- "/Users/m.hossein_drn/Documents/aim_1/final_outputs/ChromHMM_outputs/"

t_b_files_names <- c("act_cd4_rmp_hm_labels.txt", "act_cd8_rmp_hm_labels.txt", "b_mem_rmp_hm_labels.txt", "b_nai_rmp_hm_labels.txt",
                     "cd8_mem_rmp_hm_labels.txt", "cd4_nai_rmp_hm_labels.txt", "cd8_nai_rmp_hm_labels.txt", "cd4_reg_rmp_hm_labels.txt")

t_b_df_names <- c("act_cd4_t_hm", "cyto_cd8_t_hm", "mem_b_hm", "naive_b_hm", "mem_cd8_t_hm", 
                  "naive_cd4_t_hm", "naive_cd8_t_hm", "tReg_hm") # t and b cells

for (i in seq_along(t_b_files_names)) {
  file_path <- paste0(base_dir, t_b_files_names[i])
  df_name <- t_b_df_names[i]
  df <- read.delim(file_path, sep = "\t", header = TRUE)
  assign(df_name, df)
}

dc_hm = read.delim("/Users/m.hossein_drn/Documents/aim_1/final_outputs/ChromHMM_outputs/dc_rmp_hm_labels.txt", 
                   sep = "\t", header = T) # dendritic cells

mono_hm = read.delim("/Users/m.hossein_drn/Documents/aim_1/final_outputs/ChromHMM_outputs/mono_rmp_hm_labels.txt", 
                     sep = "\t", header = T) # monocytes 

nk_hm = read.delim("/Users/m.hossein_drn/Documents/aim_1/final_outputs/ChromHMM_outputs/nk_rmp_hm_labels.txt", 
                   sep = "\t", header = T) # natural killer cells 

peak_cluster_matrix = read.delim("/Users/m.hossein_drn/Documents/aim_1/final_outputs/scATAC_outputs/peak_cluster_matrix.txt", 
                                 sep = "\t", header = T)

cDC_rmps = peak_cluster_matrix[which(peak_cluster_matrix$cDC != 0),] %>% dplyr::select(cDC)
pDC_rmps = peak_cluster_matrix[which(peak_cluster_matrix$pDC != 0),] %>% dplyr::select(pDC)
iMono_rmps = peak_cluster_matrix[which(peak_cluster_matrix$iMono != 0),] %>% dplyr::select(iMono)
ncMono_rmps = peak_cluster_matrix[which(peak_cluster_matrix$ncMono != 0),] %>% dplyr::select(ncMono)
cMono_rmps = peak_cluster_matrix[which(peak_cluster_matrix$cMono != 0),] %>% dplyr::select(cMono)
adaptive_NK_rmps = peak_cluster_matrix[which(peak_cluster_matrix$adaptive_NK != 0),] %>% dplyr::select(adaptive_NK)
cyto_nk_rmps = peak_cluster_matrix[which(peak_cluster_matrix$cyto_nk != 0),] %>% dplyr::select(cyto_nk)

cDC_hm = dc_hm[which(dc_hm$rmp %in% row.names(cDC_rmps)),]
pDC_hm = dc_hm[which(dc_hm$rmp %in% row.names(pDC_rmps)),]
iMono_hm = mono_hm[which(mono_hm$rmp %in% row.names(iMono_rmps)),]
ncMono_hm = mono_hm[which(mono_hm$rmp %in% row.names(ncMono_rmps)),]
cMono_hm = mono_hm[which(mono_hm$rmp %in% row.names(cMono_rmps)),]
adaptive_NK_hm = nk_hm[which(nk_hm$rmp %in% row.names(adaptive_NK_rmps)),]
cyto_nk_hm = nk_hm[which(nk_hm$rmp %in% row.names(cyto_nk_rmps)),]

cells_without_hm = c("mkc", "plasma")

for (i in 1:nrow(final_output)) {
  cell_subtype <- final_output$cell_type[i]
  if (cell_subtype %in% cells_without_hm) {
    next
  }
  match_hm_name <- paste0(cell_subtype, "_hm")
  match_hm <- get(match_hm_name)
  rmp <- final_output$rmp[i]
  if (rmp %in% match_hm$rmp) {
    final_output$rmp_hm_label[i] <- match_hm$state[match_hm$rmp == rmp]
  }
}

# adding eQTL data
process_and_summarize_eqtl_files <- function(base_dir, file_names, df_names) {
  for (i in seq_along(file_names)) {
    file_path <- paste0(base_dir, file_names[i])
    df_name <- df_names[i]
    df <- read.delim(file_path, header = TRUE)
    df <- df %>% distinct()
    df$external_gene_name[df$external_gene_name == ""] <- df$gene_id[df$external_gene_name == ""]
    df <- df %>%
      group_by(rsid) %>%
      summarise(
        external_gene_name = paste(unique(external_gene_name), collapse = ", "),
        gene_id = paste(unique(gene_id), collapse = ", "),
        pvalue = min(pvalue)
      )
    assign(df_name, df, envir = .GlobalEnv)
  }
}

base_dir <- "/Users/m.hossein_drn/Documents/aim_1/final_outputs/eQTL_processed/"

file_names <- c("act_cd4_t_eqtl_processed_results.txt", "cd16_Mono_eqtl_processed_results.txt", 
                "cMono_eqtl_processed_results.txt", "cyto_cd8_t_eqtl_processed_results.txt",
                "mem_b_eqtl_processed_results.txt", "mem_cd8_t_eqtl_processed_results.txt", 
                "naive_b_eqtl_processed_results.txt", "naive_cd8_t_eqtl_processed_results.txt",
                "nk_eqtl_processed_results.txt", "tReg_eqtl_processed_results.txt", 
                "naive_cd4_t_eqtl_processed_results.txt")

df_names <- c("act_cd4_t_eqtl", "cd16_Mono_eqtl", "cMono_eqtl", "cyto_cd8_t_eqtl",
              "mem_b_eqtl", "mem_cd8_t_eqtl", "naive_b_eqtl", "naive_cd8_t_eqtl",
              "nk_eqtl", "tReg_eqtl", "naive_cd4_t_eqtl")

process_and_summarize_eqtl_files(base_dir, file_names, df_names)

cell_subtype_to_eqtl <- list(
  act_cd4_t = "act_cd4_t_eqtl",
  naive_cd4_t = "naive_cd4_t_eqtl",
  cyto_cd8_t = "cyto_cd8_t_eqtl",
  naive_cd8_t = "naive_cd8_t_eqtl",
  mem_cd8_t = "mem_cd8_t_eqtl",
  tReg = "tReg_eqtl",
  mem_b = "mem_b_eqtl",
  naive_b = "naive_b_eqtl",
  adaptive_NK = "nk_eqtl",
  cyto_nk = "nk_eqtl",
  ncMono = "cd16_Mono_eqtl",
  cMono = "cMono_eqtl",
  iMono = "cd16_Mono_eqtl"
)

cells_without_eqtl = c("mkc", "plasma", "cDC", "pDC")

for (i in 1:nrow(final_output)) {
  cell_subtype <- final_output$cell_type[i]
  if (cell_subtype %in% cells_without_eqtl) {
    next
  }
  match_eqtl_name <- cell_subtype_to_eqtl[[cell_subtype]]
  if (is.null(match_eqtl_name)) {
    next
  }
  match_eqtl <- get(match_eqtl_name)
  snp <- final_output$effect_snp[i]
  if (snp %in% match_eqtl$rsid) {
    final_output$eqtl_gene[i] <- match_eqtl$external_gene_name[match_eqtl$rsid == snp]
  }
}

final_output <- final_output %>%
  dplyr::select(locus, trait, loci_region, lead, lead_ppa, effect_snp, effect_ppa, tf, 
         rmp, rmp_hm_label, rmp_ppa, cicero_interacting_prom, hic_interaction_proximity, 
         directly_mapped_gene, cicero_gene, hic_gene, eqtl_gene, cell_type)

final_output[final_output == "-"] <- NA

write.csv(final_output, "complete_table.csv", row.names = FALSE)

cell_types <- unique(final_output$cell_type)
cell_type_dfs <- list()
for (cell_type in cell_types) {
  cell_type_dfs[[cell_type]] <- final_output[final_output$cell_type == cell_type, ]
}

save(cell_type_dfs, file = "cell_type_dfs.RData")

load("cell_type_dfs.RData")
View(cell_type_dfs$act_cd4_t)

