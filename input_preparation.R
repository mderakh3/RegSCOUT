setwd("/Users/m.hossein_drn/Documents/aim_1/input_preparation/")

library(xlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(GenomicRanges)
library(readr)

# read Mark et al. fine-mapping loci head
mark_fine_loci_head <- read_excel("/Users/m.hossein_drn/Documents/aim_1/input_preparation/mark_finemap.xlsx", sheet = 1)

mark_fine_loci_head <- mark_fine_loci_head[mark_fine_loci_head$signal == 1,] # remove rows with signal != 1

mark_fine_loci_head <- mark_fine_loci_head[,c("chr", "variant.lead", "position.lead", "trait.reassigned", 
                                    "credible_start", "credible_end", "all.variant")]

colnames(mark_fine_loci_head) <- c("CHR", "LeadSNP", "Position", "Trait", "credible_start", "credible_end", "credibleset")

# read the Lange et al. fine-mapping loci head
delange_fine_loci_head <- read_excel("/Users/m.hossein_drn/Documents/aim_1/input_preparation/delange_finemap.xlsx", sheet = 1)

delange_fine_loci_head <- delange_fine_loci_head[-1,]
colnames(delange_fine_loci_head) <- delange_fine_loci_head[1,]
delange_fine_loci_head <- delange_fine_loci_head[-1,]
colnames(delange_fine_loci_head)[11] <- "credibleset"

delange_fine_loci_head <- delange_fine_loci_head[,c("Chromosome", "BestSNP", "Position (bp)", "Phenotype", "credibleset")]
colnames(delange_fine_loci_head) <- c("CHR", "LeadSNP", "Position", "Trait", "credibleset")

# reading all ci-snps in delange data
files_dir <- "/Users/m.hossein_drn/Documents/aim_1/input_preparation/delange_finemap_results"
file_paths <- list.files(files_dir, pattern = "\\.finemap$", full.names = TRUE)
delange_ci <- data.frame()

delange_ci <- file_paths %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "source")

head(delange_ci)

# creating a comprehensive loci head
delange_fine_loci_head <- delange_fine_loci_head %>% separate_rows(credibleset, sep = ",")
delange_fine_loci_head <- delange_fine_loci_head[grepl("^rs", delange_fine_loci_head$credibleset),]
delange_fine_loci_head$position <- NA

for(i in 1:nrow(delange_fine_loci_head)) {
  delange_fine_loci_head$position[i] <- delange_ci$pos[delange_ci$SNP == delange_fine_loci_head$credibleset[i]]
}

delange_fine_loci_head <- delange_fine_loci_head %>% 
  mutate(credible_start = NA, credible_end = NA)

for (i in 1:nrow(delange_fine_loci_head)) {
  best_snp <- delange_fine_loci_head$LeadSNP[i]
  same_snp_rows <- delange_fine_loci_head %>% filter(LeadSNP == best_snp)
  credible_positions <- delange_fine_loci_head %>%
    filter(credibleset %in% same_snp_rows$credibleset) %>%
    pull(position)
  
  credible_start <- min(credible_positions, na.rm = TRUE)
  credible_end <- max(credible_positions, na.rm = TRUE)

  delange_fine_loci_head$credible_start[i] <- credible_start
  delange_fine_loci_head$credible_end[i] <- credible_end
}

head(delange_fine_loci_head)

delange_fine_loci_head <- delange_fine_loci_head %>%
  group_by(LeadSNP, CHR, Position, Trait, credible_start, credible_end) %>%
  summarise(
    credibleset = paste(credibleset, collapse = ", "),
    .groups = "drop"
  )

fine_loci_head <- rbind(mark_fine_loci_head, delange_fine_loci_head)

fine_loci_head <- fine_loci_head[order(fine_loci_head$CHR, fine_loci_head$credible_start),]

# defining loci and checking overlaps between ci_snps
fine_loci_head$locus <- 1:nrow(fine_loci_head)

overlaps <- c()

for (i in 1:(nrow(fine_loci_head) - 1)) {
  range1 <- GRanges(seqnames = fine_loci_head$CHR[i], 
                    ranges = IRanges(start = fine_loci_head$credible_start[i], 
                                     end = fine_loci_head$credible_end[i]))
  range2 <- GRanges(seqnames = fine_loci_head$CHR[i + 1], 
                    ranges = IRanges(start = fine_loci_head$credible_start[i + 1], 
                                     end = fine_loci_head$credible_end[i + 1]))
  if (any(overlapsAny(range1, range2))) {
    overlaps <- c(overlaps, fine_loci_head$locus[i], fine_loci_head$locus[i + 1])
  }
}

print(overlaps)

# check to see if there is shared snps in overlapping loci using credibleset and intersect function
overlapping_loci <- fine_loci_head[fine_loci_head$locus %in% overlaps,]
credibleset_56 <- unlist(strsplit(overlapping_loci$credibleset[overlapping_loci$locus == 56], ","))
credibleset_57 <- unlist(strsplit(overlapping_loci$credibleset[overlapping_loci$locus == 57], ","))

shared_snps <- intersect(credibleset_56, credibleset_57)
print(shared_snps)

library(writexl)
write_xlsx(fine_loci_head, "/Users/m.hossein_drn/Documents/aim_1/input_preparation/fine_loci_head.xlsx")

# creating ci_snps data frame for delange data
delange_ci <- delange_ci %>% filter(!is.na(P_CAUSAL))

delange_ci <- delange_ci %>% 
  mutate(A0 = NA, A1 = NA)

for (i in 1:nrow(delange_ci)) { # defining reference alleles (A0)
  if (delange_ci$beta[i] < 0) {
    delange_ci$A0[i] <- delange_ci$allele_B[i]
  } else {
    delange_ci$A0[i] <- delange_ci$allele_A[i]
  }
}

for (i in 1:nrow(delange_ci)) { # defining risk alleles (A1)
  if (delange_ci$beta[i] > 0) {
    delange_ci$A1[i] <- delange_ci$allele_B[i]
  } else {
    delange_ci$A1[i] <- delange_ci$allele_A[i]
  }
}

delange_fine_loci_head <- delange_fine_loci_head %>% separate_rows(credibleset, sep = ", ")


delange_ci_snp <- data.frame("id" = NA, "trait" = NA, "chr" = NA, "pos" = NA, 
                             "A0" = NA, "A1" = NA, "PPA" = NA)
delange_ci_snp <- delange_ci_snp[1:nrow(delange_fine_loci_head),]


for (i in 1:nrow(delange_fine_loci_head)) {
  delange_ci_snp$id[i] <- delange_fine_loci_head$credibleset[i]
  delange_ci_snp$trait[i] <- delange_fine_loci_head$Trait[i]
  delange_ci_snp$chr[i] <- delange_fine_loci_head$CHR[i]
  delange_ci_snp$pos[i] <- delange_ci$pos[delange_ci$SNP == delange_fine_loci_head$credibleset[i]]
  delange_ci_snp$A0[i] <- delange_ci$A0[delange_ci$SNP == delange_fine_loci_head$credibleset[i]]
  delange_ci_snp$A1[i] <- delange_ci$A1[delange_ci$SNP == delange_fine_loci_head$credibleset[i]]
  delange_ci_snp$PPA[i] <- delange_ci$P_CAUSAL[delange_ci$SNP == delange_fine_loci_head$credibleset[i]]
}

delange_ci_snp$author <- "delange"

# creating ci_snps data frame for mark data
mark_fine <- read_excel("/Users/m.hossein_drn/Documents/aim_1/input_preparation/mark_finemap.xlsx", sheet = 2)

mark_fine <- mark_fine[mark_fine$signal == 1,]

mark_fine <- mark_fine[,c("variant", "trait.reassigned", "chr", "position", "A0", "A1", "P_mean_95")]

mark_ci_snp <- data.frame("id" = mark_fine$variant, "trait" = mark_fine$trait.reassigned, 
                            "chr" = mark_fine$chr, "pos" = mark_fine$position, 
                            "A0" = mark_fine$A0, "A1" = mark_fine$A1, "PPA" = mark_fine$P_mean_95)

mark_ci_snp$author <- "mark"

ci_snp <- rbind(delange_ci_snp, mark_ci_snp)
rownames(ci_snp) <- NULL

ci_snp <- ci_snp[grepl("^rs", ci_snp$id),] # remove rows that their id do not start with rs

ci_snp$A0 <- toupper(ci_snp$A0) # replace every small character in A0 and A1 with the corresponding capital letter
ci_snp$A1 <- toupper(ci_snp$A1)

ci_snp$chr <- paste0("chr", ci_snp$chr) # add chr to the numbers in "chr" column

write.table(ci_snp, "/Users/m.hossein_drn/Documents/aim_1/input_preparation/ci_snp.txt", sep = "\t", quote = FALSE, row.names = FALSE)


