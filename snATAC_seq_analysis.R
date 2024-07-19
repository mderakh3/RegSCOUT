setwd("/Users/m.hossein_drn/Documents/aim_1/snATAC_seq_analysis/")

library(data.table)
library(tidyverse)
library(readxl)
library(rtracklayer)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# overlapping the peaks with the SNPs
peaks <- "/Users/m.hossein_drn/Documents/aim_1/snATAC_seq_analysis/paper/journal.pgen.1010759.s019.xlsx"
peaks_data <- read_excel(peaks)

region_ranges = IRanges(start = as.integer(peaks_data$start), 
                        end = as.integer(peaks_data$end))
region_granges = GRanges(seqnames = peaks_data$chr, ranges = region_ranges)

effect_snp <- read.delim("Ci_effect_SNPs.txt", header=TRUE, sep = "\t")

ppa_data = read.delim("ci_snp.txt", header = T, sep = "\t")

effect_snp$PPA = NA
for (i in c(1:nrow(effect_snp))){
  temp_snp = effect_snp$SNP[i]
  temp_ppa = ppa_data$PPA[ppa_data$id == temp_snp]
  if (length(temp_ppa) > 0){
    effect_snp$PPA[i] = temp_ppa
  }
}

effect_ranges = IRanges(start = effect_snp$Pos - 1, end = effect_snp$Pos)
effect_granges = GRanges(seqnames = effect_snp$CHR, ranges = effect_ranges)

reg_snp_overlap = findOverlaps(effect_granges, region_granges)
affected_regions_index = unique(subjectHits(reg_snp_overlap))

region_index = subjectHits(reg_snp_overlap)
snp_index = queryHits(reg_snp_overlap)

affected_regions = peaks_data[affected_regions_index,]

pbmc_regions = paste0(peaks_data$chr,"-", peaks_data$start, "-", peaks_data$end)

region_ppa_snp = data.frame(matrix(data = NA, nrow = length(region_index), ncol = 4)) # adding PPA information
colnames(region_ppa_snp) = c("region", "SNP", "PPA", "cell")
for (i in c(1:length(region_index))){
  temp_region = pbmc_regions[region_index[i]]
  temp_cell = peaks_data$`cell sub-types`[region_index[i]]
  temp_snp  = effect_snp$SNP[snp_index[i]]
  temp_ppa = effect_snp$PPA[snp_index[i]]
  region_ppa_snp[i,]$region = temp_region
  region_ppa_snp[i,]$cell = temp_cell
  region_ppa_snp[i,]$SNP = temp_snp
  region_ppa_snp[i,]$PPA = temp_ppa
}

region_ppa_snp <- region_ppa_snp %>% distinct(region, SNP, .keep_all = TRUE)

region_ppa_snp$sumPPA = 0 # adding the sumPPA 
for (i in c(1:nrow(region_ppa_snp))){
  temp_region = region_ppa_snp$region[i]
  temp_sumPPA = sum(region_ppa_snp$PPA[region_ppa_snp$region == temp_region])
  region_ppa_snp$sumPPA[i] = temp_sumPPA
}

region_ppa_snp <- region_ppa_snp %>% dplyr::select(region, SNP, sumPPA, cell) %>% distinct()

write.xlsx(region_ppa_snp, "risk_regions_ppa.xlsx", rowNames = FALSE)

# mapping TFs to the peaks
log_lik_ratio = read.table("log_lik_ratio.txt", sep = "\t", header = TRUE)
log_lik_ratio = log_lik_ratio %>% dplyr::select(log_lik_ratio) %>%
  mutate(log_lik_ratio = abs(log_lik_ratio)) %>%
  as.matrix()

region_tf_ratio = data.frame(matrix(data = NA, nrow = length(region_index), ncol = 5))
colnames(region_tf_ratio) = c("region", "SNP", "TF", "Pval", "cell")
for (i in c(1:length(region_index))){
  temp_region = pbmc_regions[region_index[i]]
  temp_cell = peaks_data$`cell sub-types`[region_index[i]]
  temp_tf = effect_snp$TF[snp_index[i]]
  temp_snp  = effect_snp$SNP[snp_index[i]]
  temp_pval = effect_snp$FDR_corr_P_value[snp_index[i]]
  region_tf_ratio[i,]$region = temp_region
  region_tf_ratio[i,]$cell = temp_cell
  region_tf_ratio[i,]$TF = temp_tf
  region_tf_ratio[i,]$SNP = temp_snp
  region_tf_ratio[i,]$Pval = temp_pval
}

region_tf_ratio$TFSNP = paste0(region_tf_ratio$TF, "-", region_tf_ratio$SNP)
region_tf_ratio$log_lik_ratio = 0
for (i in c(1:nrow(region_tf_ratio))){
  temp_TFSNP = region_tf_ratio$TFSNP[i]
  temp_log_lik_ratio = log_lik_ratio[temp_TFSNP,]
  region_tf_ratio$log_lik_ratio[i] = temp_log_lik_ratio
} 
region_tf_ratio <- region_tf_ratio %>% dplyr::select(region, cell, TFSNP, log_lik_ratio) %>% distinct()

write.xlsx(region_tf_ratio, "risk_regions_ratio.xlsx", rowNames = FALSE)

# rmp heatmap and matrix
output_file_main = "/Users/m.hossein_drn/Documents/aim_1/snATAC_seq_analysis/"

rmp_matrix <- region_ppa_snp %>%
  dplyr::select(region, cell) %>% distinct() %>%
  separate_rows(cell, sep = ",")

rmp_matrix <- rmp_matrix %>%
  mutate(value = 1) %>%
  spread(key = cell, value = value, fill = 0) %>% 
  tibble::column_to_rownames(var = "region") %>%
  as.matrix()

ppa_matrix <- region_ppa_snp %>%
  dplyr::select(region, sumPPA) %>% 
  distinct() %>%
  tibble::column_to_rownames(var = "region") %>%
  as.matrix()

ppa_matrix <- ppa_matrix[rownames(rmp_matrix), , drop = FALSE]

f1 = colorRamp2(seq(0, 1, length = 2), c("#EEEEEE", "red"))
heatmap_peaks <- Heatmap(rmp_matrix, name = "Percentage accessibility", col = f1, 
                         column_title = "Immune Cell subtypes", row_title = "Risk-mediating Peaks",
                         row_names_gp = grid::gpar(fontsize = 15),
                         column_names_gp = grid::gpar(fontsize = 15),
                         rect_gp = gpar(col= "#84878a"),
                         heatmap_legend_param = list(title = "Accessibility", at = c(0, 1), 
                                                     labels = c("0", "1"), 
                                                     color_bar = "vertical", 
                                                     legend_width = unit(16, "cm")))

f2 = colorRamp2(seq(0, max(ppa_matrix[,1], na.rm = TRUE), length = 2), c("#EEEEEE","blue"))
heatmap_ppa <- Heatmap(ppa_matrix, name = "SumPPA", col = f2, 
                       heatmap_legend_param = list(title = "SumPPA", 
                                                   at = c(0, max(ppa_matrix, na.rm = TRUE)), 
                                                   labels = c("0", round(max(ppa_matrix, na.rm = TRUE))), 
                                                   color_bar = "vertical", 
                                                   legend_width = unit(16, "cm")),
                       row_names_gp = grid::gpar(fontsize = 4),
                       column_names_gp = grid::gpar(fontsize = 15),
                       layer_fun = function(j, i, x, y, width, height, fill) {
                         grid::grid.text(sprintf("%.3f", ppa_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 4))
                       })

output_peak_file = paste0(output_file_main, "cell_peak.png")
png(output_peak_file, width = 2400, height = 8000, res = 300)
combined_heatmaps <- HeatmapList(heatmap_peaks + heatmap_ppa)
print(combined_heatmaps)
dev.off()

write.table(rmp_matrix, "peak_cluster_matrix.txt", sep = "\t", row.names = TRUE)

# tf heatmap and matrix
output_tf_file = paste0(output_file_main, "cell_tf.png")

tf_matrix <- region_tf_ratio %>%
  dplyr::select(TFSNP, cell) %>% distinct() %>%
  separate_rows(cell, sep = ",") %>%
  mutate(value = 1) %>%
  spread(key = cell, value = value, fill = 0) %>%
  tibble::column_to_rownames(var = "TFSNP") %>%
  as.matrix()

llr_matrix <- region_tf_ratio %>%
  dplyr::select(TFSNP, log_lik_ratio) %>% distinct() %>%
  tibble::column_to_rownames(var = "TFSNP") %>%
  as.matrix()

llr_matrix <- llr_matrix[rownames(tf_matrix), , drop = FALSE]

f1 = colorRamp2(seq(0, 1, length = 2), c("#EEEEEE", "red"))
heatmap_tf <- Heatmap(tf_matrix, name = "Pair Presence", col = f1, 
                      column_title = "Immune Cell Subtypes", row_title = "TF-SNP",
                      row_names_gp = grid::gpar(fontsize = 15),
                      column_names_gp = grid::gpar(fontsize = 15),
                      rect_gp = gpar(col= "#84878a"),
                      heatmap_legend_param = list(title = "TF Presence", at = c(0, 1), 
                                                  labels = c("0", "1"), 
                                                  color_bar = "vertical", 
                                                  legend_width = unit(16, "cm")))

f2 = colorRamp2(seq(0, max(llr_matrix[,1], na.rm = TRUE), length = 2), c("#EEEEEE","green"))
heatmap_llr <- Heatmap(llr_matrix, name = "LLR", col = f2, 
                       heatmap_legend_param = list(title = "abs(LLR)", 
                                                   at = c(0, max(llr_matrix, na.rm = TRUE)), 
                                                   labels = c("0", round(max(llr_matrix, na.rm = TRUE))), 
                                                   color_bar = "vertical", 
                                                   legend_width = unit(16, "cm")),
                       row_names_gp = grid::gpar(fontsize = 2),
                       column_names_gp = grid::gpar(fontsize = 15),
                       layer_fun = function(j, i, x, y, width, height, fill) {
                         grid::grid.text(sprintf("%.3f", llr_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 2))
                       })

png(output_tf_file, width = 2400, height = 8000, res = 300)
combined_heatmaps <- HeatmapList(heatmap_tf + heatmap_llr)
print(combined_heatmaps)
dev.off()

write.table(tf_matrix, "tf_cluster_matrix.txt", sep = "\t", row.names = TRUE)

# finding proximal RMPs
library(ape)
gene_annot_dir = "/Users/m.hossein_drn/Documents/aim_1/snATAC_seq_analysis/gencode.v19.annotation.gff3"
gene_annot = read.gff(gene_annot_dir, na.strings = c(".", "?"), GFF3 = TRUE)

gene_type_list = str_split(gene_annot$attributes, "gene_type=", simplify = TRUE)
gene_type_list = str_split(gene_type_list[,2], ";", simplify = TRUE)[,1]

gene_coding_index = gene_type_list == "protein_coding"
gene_annot = gene_annot[gene_coding_index,]

gene_transcript_data = gene_annot[gene_annot$type == "transcript",]

gene_id_list = gene_transcript_data$attributes
head(gene_id_list)

genome_built = "hg19"

if (genome_built == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome = "BSgenome.Hsapiens.UCSC.hg19"  
  hg_19_ref = BSgenome.Hsapiens.UCSC.hg19@seqinfo
  hg_19_dat = as.data.frame(matrix(NA, nrow = length(hg_19_ref@seqnames), ncol = 2))
  hg_19_dat$V1 = hg_19_ref@seqnames
  hg_19_dat$V2 = hg_19_ref@seqlengths
}else if(genome_built == "hg38"){
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome = "BSgenome.Hsapiens.UCSC.hg38" 
  hg_38_ref = BSgenome.Hsapiens.UCSC.hg38@seqinfo
  hg_38_dat = as.data.frame(matrix(NA, nrow = length(hg_38_ref@seqnames), ncol = 2))
  hg_38_dat$V1 = hg_38_ref@seqnames
  hg_38_dat$V2 = hg_38_ref@seqlengths
}

gene_id_list = str_split(gene_id_list, "gene_name=", simplify = TRUE)
head(gene_id_list,1)

gene_id_list = str_split(gene_id_list[,2], ";", simplify = TRUE)[,1]
print(head(gene_id_list,10))

gene_transcript_data[["gene_name"]] = gene_id_list

pos_strand_index = gene_transcript_data$strand == "+"
neg_strand_index = gene_transcript_data$strand == "-"

gene_transcript_data[["TSS"]] = rep(0, nrow(gene_transcript_data))
gene_transcript_data$TSS[pos_strand_index] = gene_transcript_data$start[pos_strand_index]
gene_transcript_data$TSS[neg_strand_index] = gene_transcript_data$end[neg_strand_index]

pos_gene_ranges = IRanges(start = gene_transcript_data$TSS[pos_strand_index] - 2000, 
                          end = gene_transcript_data$TSS[pos_strand_index] + 2000)
pos_gene_granges = GRanges(ranges = pos_gene_ranges, seqnames = gene_transcript_data$seqid[pos_strand_index],
                           strand = gene_transcript_data$strand[pos_strand_index])
pos_gene_granges$gene_name = gene_transcript_data$gene_name[pos_strand_index]

neg_gene_ranges = IRanges(start = gene_transcript_data$TSS[neg_strand_index] - 2000, 
                          end = gene_transcript_data$TSS[neg_strand_index] + 2000)
neg_gene_granges = GRanges(ranges = neg_gene_ranges, seqnames = gene_transcript_data$seqid[neg_strand_index],
                           strand = gene_transcript_data$strand[neg_strand_index])
neg_gene_granges$gene_name = gene_transcript_data$gene_name[neg_strand_index]

gene_granges = c(pos_gene_granges, neg_gene_granges)

pbmc_promoters = gene_granges

print("Gene data loaded successfully")

rmp_ranges = IRanges(start = as.integer(affected_regions$start), end = as.integer(affected_regions$end))
rmp_granges = GRanges(seqnames = affected_regions$chr, ranges = rmp_ranges)

overlap = findOverlaps(rmp_granges, pbmc_promoters)
pbmc_promoters[subjectHits(overlap)]
rmp_granges[queryHits(overlap)]

rmp_prom_data = data.frame(rmp_granges[queryHits(overlap)], pbmc_promoters[subjectHits(overlap)])
colnames(rmp_prom_data) = c("rmp_chr", "rmp_start", "rmp_end", "rmp_width", "rmp_strand",
                            "prom_chr", "prom_start", "prom_end", "prom_width", "prom_strand",
                            "gene_name")

rmp_prom_data$rmp_prom = paste0(rmp_prom_data$rmp_chr,"-", rmp_prom_data$rmp_start, "-", rmp_prom_data$rmp_end)
rmp_prom_data = rmp_prom_data %>% 
  select(rmp_prom, gene_name, rmp_strand, 
         prom_strand, rmp_width, prom_width) %>% 
  distinct()

write.xlsx(rmp_prom_data, file = "rmp_promoter_data.xlsx", sheetName = "rmp_promoter_data")

# finding distal RMPs
cicero_matrix <- read.table("annotated.coaccessible.matrix4s", header = TRUE, sep = "\t")
cicero_matrix <- cicero_matrix %>% select(-loops)

cicero_matrix$Peak1 <- gsub("_", "-", cicero_matrix$Peak1)
cicero_matrix$Peak2 <- gsub("_", "-", cicero_matrix$Peak2)

for (i in 1:nrow(cicero_matrix)) { # defining promoter column
  if (cicero_matrix$prom_side[i] == "A") {
    cicero_matrix$promoter[i] <- cicero_matrix$Peak1[i]
  } else {
    cicero_matrix$promoter[i] <- cicero_matrix$Peak2[i]
  }
}

for (i in 1:nrow(cicero_matrix)) { # defining enhancer column
  if (cicero_matrix$prom_side[i] == "B") {
    cicero_matrix$enhancer[i] <- cicero_matrix$Peak1[i]
  } else {
    cicero_matrix$enhancer[i] <- cicero_matrix$Peak2[i]
  }
}

cic_rmp = cicero_matrix %>% filter(enhancer %in% rownames(rmp_matrix))

write.xlsx(cic_rmp, "cic_rmp_enhancers.xlsx", rowNames = FALSE)
