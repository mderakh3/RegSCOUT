setwd("/Users/m.hossein_drn/Documents/aim_1/histone_marks_analysis/")

library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(stringr)

# preparing ChromHMM annotation files
main_dir = "/Users/m.hossein_drn/Documents/aim_1/histone_marks_analysis/ChromHMM" 
chromhmm_datasets <- read.delim("chromhmm_datasets.txt", header = TRUE, sep = "\t")
chain_dir = "/Users/m.hossein_drn/Documents/aim_1/histone_marks_analysis/hg38ToHg19.over.chain"
chain_data = import.chain(chain_dir)

sub_dirs <- c("dendritic_cells", "monocytes", "b_mem_cells", "b_nai_cells", "act_cd4_cells", 
              "act_cd8_cells", "cd4_nai_cells", "cd4_reg_cells", 
              "cd8_nai_cells", "cd8_mem_cells", "nk_cells")

file_names <- c("ENCFF022ZIV.bed", "ENCFF269WBG.bed", "ENCFF476DBW.bed", "ENCFF478AZW.bed", 
                "ENCFF342HTO.bed", "ENCFF555HRT.bed", "ENCFF462ZAZ.bed", "ENCFF255WPS.bed", 
                "ENCFF890OEH.bed", "ENCFF358UZU.bed", "ENCFF953PSB.bed")

cell_types <- c("dc_his", "mono_his", "b_mem_his", "b_nai_his", "act_cd4_his", 
                "act_cd8_his", "cd4_nai_his", "cd4_reg_his", 
                "cd8_nai_his", "cd8_mem_his", "nk_his")

his_dirs <- setNames(vector("list", length(cell_types)), cell_types)

for (i in seq_along(sub_dirs)) {
  his_dirs[[cell_types[i]]] <- paste0(main_dir, "/", sub_dirs[i], "/", file_names[i])
}

his_data_list <- setNames(vector("list", length(cell_types)), cell_types)

col_names <- c("chr", "start", "end", "state", "score", "strand", "thickStart", "thickEnd", "RGB")

for (i in seq_along(his_dirs)) {
  his_data_list[[i]] <- fread(his_dirs[[i]], header = FALSE, sep = "\t")
  setnames(his_data_list[[i]], col_names)
}

for (i in seq_along(his_data_list)) {
  cat(paste("Data for", names(his_data_list)[i], ":\n"))
  print(head(his_data_list[[i]]))
  cat("\n")
}

# creating GRanges objects for histone marks data
his_granges <- list()

for (cell in names(his_data_list)) {
  if (chromhmm_datasets$genome_built[chromhmm_datasets$cell_type == cell] == "hg19") {
    his_granges[[cell]] <- GRanges(seqnames = his_data_list[[cell]]$chr, ranges = IRanges(start = his_data_list[[cell]]$start, end = his_data_list[[cell]]$end))
    his_granges[[cell]]$state = his_data_list[[cell]]$state
  } else if (chromhmm_datasets$genome_built[chromhmm_datasets$cell_type == cell] == "hg38") {
    his_granges[[cell]] <- GRanges(seqnames = his_data_list[[cell]]$chr, ranges = IRanges(start = his_data_list[[cell]]$start, end = his_data_list[[cell]]$end))
    his_granges[[cell]]$state = his_data_list[[cell]]$state
    his_granges[[cell]] <- unlist(liftOver(his_granges[[cell]], chain_data))
  }
}

peak_cluster_matrix = read.delim("peak_cluster_matrix.txt", header = TRUE, sep = "\t") # preparing RMP data

cell_types <- list(
  dc_rmp = c("cDC", "pDC"),
  mono_rmp = c("cMono", "ncMono", "iMono"),
  b_mem_rmp = "mem_b",
  b_nai_rmp = "naive_b",
  act_cd4_rmp = "act_cd4_t",
  act_cd8_rmp = "cyto_cd8_t", 
  cd4_nai_rmp = "naive_cd4_t",
  cd4_reg_rmp = "tReg",
  cd8_nai_rmp = "naive_cd8_t",
  cd8_mem_rmp = "mem_cd8_t",
  nk_rmp = c("adaptive_NK", "cyto_nk")
) 

# preparing cell subtype specific RMP data
result_matrices <- list()

for (cell_type in names(cell_types)) {
  relevant_columns <- cell_types[[cell_type]]
  sub_matrix <- peak_cluster_matrix[, relevant_columns, drop = FALSE]
  filtered_matrix <- sub_matrix[rowSums(sub_matrix) != 0, , drop = FALSE]
  result_matrices[[cell_type]] <- filtered_matrix
}

for (i in seq_along(result_matrices)) {
  result_matrices[[i]]$rmp <- rownames(result_matrices[[i]])
  rownames(result_matrices[[i]]) <- NULL
} 

# creating GRanges objects for RMP data
rmp_granges <- list()

for (cell_type in names(result_matrices)) {
  rmp_granges[[cell_type]] <- GRanges(seqnames = sapply(result_matrices[[cell_type]]$rmp, function(x) str_split(x, "-")[[1]][1]), 
                                      ranges = IRanges(start = as.numeric(sapply(result_matrices[[cell_type]]$rmp, function(x) str_split(x, "-")[[1]][2])), 
                                                       end = as.numeric(sapply(result_matrices[[cell_type]]$rmp, function(x) str_split(x, "-")[[1]][3]))))
}

match_cell_types <- list(
  dc_rmp = "dc_his",
  mono_rmp = "mono_his",
  b_mem_rmp = "b_mem_his",
  b_nai_rmp = "b_nai_his",
  act_cd4_rmp = "act_cd4_his",
  act_cd8_rmp = "act_cd8_his",
  cd4_nai_rmp = "cd4_nai_his",
  cd4_reg_rmp = "cd4_reg_his",
  cd8_nai_rmp = "cd8_nai_his",
  cd8_mem_rmp = "cd8_mem_his",
  nk_rmp = "nk_his"
)

# adding ChromHMM states to RMP data
rmp_labels_list <- list()

for (cell_type in names(match_cell_types)) {
  rmp_grange <- rmp_granges[[cell_type]]
  his_grange <- his_granges[[match_cell_types[[cell_type]]]]
  
  rmp_his_index <- findOverlaps(rmp_grange, his_grange)
  
  query_hits <- queryHits(rmp_his_index)
  subject_hits <- subjectHits(rmp_his_index)
  
  query_granges <- rmp_grange[query_hits]
  subject_granges <- his_grange[subject_hits]
  state_values <- mcols(subject_granges)$state
  
  rmp_labels <- data.frame(
    Query = as.character(query_granges),
    Subject = as.character(subject_granges),
    State = state_values
  )
  
  rmp_labels$Query <- gsub(":", "-", rmp_labels$Query)
  rmp_labels$Subject <- gsub(":", "-", rmp_labels$Subject)
  colnames(rmp_labels) <- c("rmp", "his", "state")
  
  rmp_labels <- rmp_labels[order(rmp_labels$rmp), ]
  rmp_labels$state <- as.character(rmp_labels$state)
  
  for (i in 1:nrow(rmp_labels)) {
    if (rmp_labels$rmp[i] %in% rmp_labels$rmp[-i]) {
      rmp_labels$state[i] <- paste(rmp_labels$state[rmp_labels$rmp == rmp_labels$rmp[i]], collapse = ", ")
      rmp_labels <- rmp_labels[-which(rmp_labels$rmp == rmp_labels$rmp[i] & rmp_labels$state != rmp_labels$state[i]), ]
    }
  }
  
  rmp_his <- merge(result_matrices[[cell_type]], rmp_labels, by = "rmp", all.x = TRUE)
  rmp_his <- rmp_his[, c("rmp", "state")]
  
  if (length(levels(as.factor(his_granges[[match_cell_types[[cell_type]]]]$state))) == 15) {
    rmp_his$version <- "ChromHMM 15-state"
  } else if (length(levels(as.factor(his_granges[[match_cell_types[[cell_type]]]]$state))) == 18) {
    rmp_his$version <- "ChromHMM 18-state"
  }
  
  rmp_labels_list[[cell_type]] <- rmp_his
}

for (cell_type in names(rmp_labels_list)) {
  write.table(rmp_labels_list[[cell_type]], file = paste0("ChromHMM_outputs/", cell_type, "_hm_labels.txt"), sep = "\t", row.names = FALSE)
}
