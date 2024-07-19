setwd("/Users/m.hossein_drn/Documents/aim_1/HiC_enhancers_evaluation/")

library(xlsx)
library(GenomicRanges)
library(stringr)
library(liftOver)
library(rtracklayer)
library(dplyr)

# preparing input HiC data
pchic_matrix_dir = "/Users/m.hossein_drn/Documents/aim_1/HiC_enhancers_evaluation/paper/DATA_S1/PCHiC_peak_matrix_cutoff5.tsv"
pchic_matrix = read.delim(file = pchic_matrix_dir, sep = "\t", header = TRUE) # pcHi-c data from Biola M. Javierre et al

pchic_matrix <- pchic_matrix %>%
  mutate(
    promoter = paste0("chr", baitChr, "-", baitStart, "-", baitEnd),
    enhancer = paste0("chr", oeChr, "-", oeStart, "-", oeEnd)
  ) %>%
  dplyr::select(promoter, enhancer, Mon, tCD4, tCD8, tB, baitName, oeName)

colnames(pchic_matrix) = c("promoter", "enhancer", "Mon", "tCD4", "tCD8", "tB", "ba_gene", "oe_gene")

inthic_loops_dir = "/Users/m.hossein_drn/Documents/aim_1/HiC_enhancers_evaluation/ENCFF492SAH.bedpe" # intHiC for NK cells
inthic_loops = read.delim(inthic_loops_dir, header = TRUE, sep = "\t")

inthic_loops <- inthic_loops %>%
  filter(grepl("^chr", X.chr1)) %>%
  filter(fdrBL < 0.05, fdrDonut < 0.05, fdrH < 0.05, fdrV < 0.05)

colnames(inthic_loops) = c("chr1", "start1", "end1", "chr2", "start2", "end2", "score")
inthic_loops = inthic_loops %>% select(chr1, start1, end1, chr2, start2, end2)

chain_dir = "/Users/m.hossein_drn/Documents/aim_1/HiC_enhancers_evaluation/hg38ToHg19.over.chain"
chain = import.chain(chain_dir)

lifted_regions <- data.frame()

# doing liftover and skipping pairs for which the liftover returns more than one interval
for (i in 1:nrow(inthic_loops)) { 
  row <- inthic_loops[i,]
  
  pair1_grange <- GRanges(seqnames = row$chr1, ranges = IRanges(start = row$start1, end = row$end1))
  pair1_lifted <- liftOver(pair1_grange, chain)
  
  pair2_grange <- GRanges(seqnames = row$chr2, ranges = IRanges(start = row$start2, end = row$end2))
  pair2_lifted <- liftOver(pair2_grange, chain)
  
  # check if both liftOvers were successful and returned exactly one interval
  if (!is.null(pair1_lifted[[1]]) && length(pair1_lifted[[1]]) == 1 &&
      !is.null(pair2_lifted[[1]]) && length(pair2_lifted[[1]]) == 1) {
    
    lifted_pair1 <- as.data.frame(pair1_lifted[[1]]) # create data frames from the lifted GRanges objects
    lifted_pair2 <- as.data.frame(pair2_lifted[[1]])
    
    loops <- cbind(lifted_pair1, lifted_pair2) # combine the lifted pairs
    lifted_regions <- rbind(lifted_regions, loops)
  }
}

colnames(lifted_regions) = c("chr1", "start1", "end1", "width1", "strand1", "chr2", "start2", "end2", "width2", "strand2")
write.table(lifted_regions, file = "hg19_nk_intHiC.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cic_pir_dir = "/Users/m.hossein_drn/Documents/aim_1/HiC_enhancers_evaluation/cic_rmp_enhancers.xlsx"
cic_pir = read.xlsx(cic_pir_dir, 1, header = TRUE, row.names = NULL)

# preparing cicero data for t, b, and mono
process_cicero <- function(data, cell_type) { 
  data %>%
    dplyr::select(promoter, enhancer, all_of(cell_type), genes) %>%
    filter(!is.na(.data[[cell_type]]))
}

cic_pir_b <- process_cicero(cic_pir, "b")
cic_pir_mono <- process_cicero(cic_pir, "mono")
cic_pir_t <- process_cicero(cic_pir, "t")

# preparing cicero data for nk
cic_pir_nk = cic_pir %>% dplyr::select(Peak1, Peak2, promoter, enhancer, nk, genes) %>%
  filter(!is.na(nk))

# preparing pcHi-C data for B, T, and Mono
process_pchic <- function(data, cell_type, new_col_name) { 
  data %>%
    dplyr::select(promoter, enhancer, all_of(cell_type), ba_gene, oe_gene) %>%
    filter(.data[[cell_type]] >= 5)
}

pchic_pir_b <- process_pchic(pchic_matrix, "tB", "tB")
pchic_pir_mono <- process_pchic(pchic_matrix, "Mon", "Mon")

pchic_pir_tCD4 <- process_pchic(pchic_matrix, "tCD4", "tT")
pchic_pir_tCD8 <- process_pchic(pchic_matrix, "tCD8", "tT")
pchic_pir_t <- bind_rows(pchic_pir_tCD4, pchic_pir_tCD8)

# creating granges for t, b, and mono
create_granges <- function(data, cell_type_prefix) { 
  prom_df <- as.data.frame(str_split(data$promoter, "-", simplify = TRUE))
  colnames(prom_df) <- c("chr", "start", "end")
  
  enh_df <- as.data.frame(str_split(data$enhancer, "-", simplify = TRUE))
  colnames(enh_df) <- c("chr", "start", "end")
  
  loops <- list(prom_granges = GRanges(seqnames = prom_df$chr, ranges = IRanges(start = as.numeric(prom_df$start), end = as.numeric(prom_df$end))),
                enh_granges = GRanges(seqnames = enh_df$chr, ranges = IRanges(start = as.numeric(enh_df$start), end = as.numeric(enh_df$end)))
  )
  
  return(loops)
}

cic_b_loops <- create_granges(cic_pir_b, "b") # B cells
pchic_b_loops <- create_granges(pchic_pir_b, "tB")

cic_mono_loops <- create_granges(cic_pir_mono, "mono") # monocytes
pchic_mono_loops <- create_granges(pchic_pir_mono, "Mon")

cic_t_loops <- create_granges(cic_pir_t, "t") # T cell
pchic_t_loops <- create_granges(pchic_pir_t, "tT")


cic_nk_peak1 = as.data.frame(str_split(cic_pir_nk$Peak1, "-", simplify = TRUE)) # NK cells
colnames(cic_nk_peak1) = c("chr", "start", "end")
cic_nk_peak2 = as.data.frame(str_split(cic_pir_nk$Peak2, "-", simplify = TRUE))
colnames(cic_nk_peak2) = c("chr", "start", "end")

cic_nk_loops = list(GRanges(seqnames = cic_nk_peak1$chr, ranges = IRanges(start = as.numeric(cic_nk_peak1$start), end = as.numeric(cic_nk_peak1$end))),
                    GRanges(seqnames = cic_nk_peak2$chr, ranges = IRanges(start = as.numeric(cic_nk_peak2$start), end = as.numeric(cic_nk_peak2$end))))

# function to overlap loops
interact_peak_overlap = function(peak_list1, peak_list2){ 
  
  list1_grgr1 = peak_list1[[1]]
  list1_grgr2 = peak_list1[[2]]
  
  list2_grgr1 = peak_list2[[1]]
  list2_grgr2 = peak_list2[[2]]
  
  grg1_overlaps = findOverlaps(list2_grgr1, list1_grgr1)
  grg2_overlaps = findOverlaps(list2_grgr2, list1_grgr2)
  
  index1 = queryHits(grg1_overlaps[grg1_overlaps %in% grg2_overlaps])
  index2 = subjectHits(grg1_overlaps[grg1_overlaps %in% grg2_overlaps])
  
  return(paste0(index1, "-", index2))
}

cell_types <- list( # for B, T, and Mono
  list(name = "b", cic_loops = cic_b_loops, pchic_loops = pchic_b_loops, cic_data = cic_pir_b, pchic_data = pchic_pir_b),
  list(name = "mono", cic_loops = cic_mono_loops, pchic_loops = pchic_mono_loops, cic_data = cic_pir_mono, pchic_data = pchic_pir_mono),
  list(name = "t", cic_loops = cic_t_loops, pchic_loops = pchic_t_loops, cic_data = cic_pir_t, pchic_data = pchic_pir_t)
)

confirmed_regions_list <- list()

main_out_dir <- "/Users/m.hossein_drn/Documents/aim_1/HiC_enhancers_evaluation/"

# Cicero overlap with pcHiC for B, T, and Mono
for (cell in cell_types) { 
  cell_name <- cell$name
  cic_loops <- cell$cic_loops
  pchic_loops <- cell$pchic_loops
  cic_data <- cell$cic_data
  pchic_data <- cell$pchic_data
  
  confirmed_regions <- interact_peak_overlap(pchic_loops, cic_loops)
  confirmed_regions <- as.data.frame(str_split(confirmed_regions, "-", simplify = TRUE))
  colnames(confirmed_regions) <- c("cicero_index", "pchic_index")
  
  cic_confirmed <- cic_data[as.numeric(confirmed_regions$cicero_index),]
  pchic_confirmed <- pchic_data[as.numeric(confirmed_regions$pchic_index),]
  
  final_data <- cbind(cic_confirmed, pchic_confirmed)
  colnames(final_data) <- c("cic_promoter", "cic_enhancer", paste0("cic_", cell_name), "cic_genes",
                            "pchic_promoter", "pchic_enhancer", paste0("pchic_", cell_name), "pchic_ba_genes", "pchic_oe_genes")
  
  write.table(final_data, file = paste0(main_out_dir, cell_name, "_pcHiC_Cic.txt"), row.names = FALSE, col.names = TRUE, sep = "\t")
  # print a message to confirm processing of the cell type
  print(paste("Processed and saved data for cell type:", cell_name))
}

# Cicero overlap with intHiC for NK
inthic_nk_loops = list(GRanges(seqnames = lifted_regions$chr1, ranges = IRanges(start = as.numeric(lifted_regions$start1), end = as.numeric(lifted_regions$end1))),
                       GRanges(seqnames = lifted_regions$chr2, ranges = IRanges(start = as.numeric(lifted_regions$start2), end = as.numeric(lifted_regions$end2))))

confirmed_regions = interact_peak_overlap(inthic_nk_loops, cic_nk_loops)

confirmed_regions = as.data.frame(str_split(confirmed_regions, "-", simplify = TRUE))
colnames(confirmed_regions) = c("cicero_index", "inthic_index")

cic_nk_confirmed = cic_pir_nk[as.numeric(confirmed_regions$cicero_index),]
colnames(cic_nk_confirmed) = c("cic_peak1", "cic_peak2", "cic_promoter", "cic_enhancer", "cic_nk", "cic_genes")

inthic_nk_confirmed = lifted_regions[as.numeric(confirmed_regions$inthic_index),]
inthic_nk_confirmed$inthic_peak1 = paste0(inthic_nk_confirmed$chr1, "-", inthic_nk_confirmed$start1, "-", inthic_nk_confirmed$end1)
inthic_nk_confirmed$inthic_peak2 = paste0(inthic_nk_confirmed$chr2, "-", inthic_nk_confirmed$start2, "-", inthic_nk_confirmed$end2)
inthic_nk_confirmed = inthic_nk_confirmed[,c("inthic_peak1", "inthic_peak2", "width1", "width2")]

nk_cells = cbind(cic_nk_confirmed, inthic_nk_confirmed)
write.table(nk_cells, file = paste0(main_out_dir, "nk_intHiC_Cic.txt"), row.names = TRUE, col.names = TRUE, sep = "\t")
