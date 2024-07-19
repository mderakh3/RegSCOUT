setwd("/Users/m.hossein_drn/Documents/aim_1/tf_binding_analysis/")

library(atSNP)
library(TFBSTools)
library(R.utils)
library(tidyr)
library(dplyr)
library(stringr)

# preparing the JASPAR PWMs
jaspar_pwm <- readJASPARMatrix("/Users/m.hossein_drn/Documents/aim_1/tf_binding_analysis/combined_755_motifs_jaspar2024.txt", matrixClass = c("PFM", "PWM", "PWMProb"))

jaspar_pwm_list = jaspar_pwm@listData

jaspar_pwm_atsnp = list()

for (i in c(1:length(jaspar_pwm_list))){
  temp_tf = jaspar_pwm_list[[i]]
  temp_tf_name = temp_tf@name
  temp_tf_matrix = temp_tf@profileMatrix
  new_matrix = apply(temp_tf_matrix,2,function(x)(x/sum(x)))
  jaspar_pwm_atsnp[[temp_tf_name]] = t(new_matrix)
}

ci_snp_dir = "/Users/m.hossein_drn/Documents/aim_1/tf_binding_analysis/ci_snp.txt"
ci_snp = read.table(ci_snp_dir, sep = "\t", header = TRUE)

# separating those rows that have more than one letter in either A0 or A1 column  
trash_ci_snp = ci_snp[which(nchar(as.character(ci_snp$A0)) != 1 | nchar(as.character(ci_snp$A1)) != 1),]
ci_snp = ci_snp[which(nchar(as.character(ci_snp$A0)) == 1 & nchar(as.character(ci_snp$A1)) == 1),]
print(paste0("Precentage of SNPs with more than one letter in A0 or A1 column: ", nrow(trash_ci_snp)/nrow(ci_snp)*100))

genome_built = "hg19"

snp_table = as.data.frame(matrix(0, nrow = nrow(ci_snp), ncol = 5))
colnames(snp_table) = c("chr","snp","snpid","a1","a2")

for (i in c(1:nrow(ci_snp))){
  temp_chr = ci_snp$chr[i]
  temp_pos = ci_snp$pos[i]
  temp_snpid = ci_snp$id[i]
  temp_a1 = ci_snp$A0[i]
  temp_a2 = ci_snp$A1[i]
  
  snp_table[i,] = c(temp_chr, temp_pos, temp_snpid, temp_a1, temp_a2)
}

snp_table_dir = "/Users/m.hossein_drn/Documents/aim_1/tf_binding_analysis/snp_table.txt"
write.table(snp_table, file = snp_table_dir, col.names = TRUE,
            row.names = FALSE, quote = FALSE)

library(BSgenome.Hsapiens.UCSC.hg19)

reg1_snp_data = LoadSNPData(filename = snp_table_dir,
                            genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                            half.window.size = 30, default.par = TRUE,
                            mutation = FALSE)

results = ComputeMotifScore(jaspar_pwm_atsnp, reg1_snp_data, ncores = 2)

all_results <- list()

for (i in 1:10) {
  results_pval_i <- ComputePValues(jaspar_pwm_atsnp, reg1_snp_data, results$motif.scores, ncores = 2, testing.mc=TRUE)
  all_results[[i]] <- results_pval_i
}

results_pval <- do.call(rbind, all_results)

results_pval_val = results_pval$pval_diff
results_pval_val_cor = p.adjust(results_pval_val, method = "fdr")
results_pval$results_pval_val_cor <- results_pval_val_cor

# counting the number of significant snp-tf pairs in all runs
data_with_counts <- results_pval %>%
  group_by(snpid, motif) %>%
  mutate(
    Count = n(),  # total count of each SNP-TF pair
    Count_FDR_lt_0.05 = sum(results_pval_val_cor < 0.05)  # count of FDR values < 0.05
  )

write.table(data_with_counts, "data_with_counts.txt", sep = "\t")

results_pval_sig = data_with_counts[data_with_counts$Count_FDR_lt_0.05 >= 2,]
results_pval_sig = results_pval_sig[results_pval_sig$results_pval_val_cor < 0.05,] 

reg1_effect_snps = sort(results_pval_sig$snpid)

final_effect_snp_frame = data.frame(matrix(data = NA, 
                                           nrow = length(results_pval_sig$snpid),
                                           ncol = 7))

colnames(final_effect_snp_frame) = c("SNP","CHR","Pos","TF", "Motif_and_TF", "Raw_P_value", "FDR_corr_P_value")

for (i in c(1:length(results_pval_sig$snpid))){
  temp_id = results_pval_sig$snpid[i]
  ci_effect_index = which(ci_snp$id == temp_id)
  temp_pos = ci_snp$pos[ci_effect_index]
  temp_chr = ci_snp$chr[ci_effect_index]
  temp_TF = substring(results_pval_sig$motif[i], 10)
  temp_TF_motif = results_pval_sig$motif[i]
  temp_pval = results_pval_sig$pval_diff[i]
  temp_pval_corr = results_pval_sig$results_pval_val_cor[i]
  
  final_effect_snp_frame[i,]$SNP = temp_id
  final_effect_snp_frame[i,]$CHR = temp_chr
  final_effect_snp_frame[i,]$Pos = temp_pos
  final_effect_snp_frame[i,]$TF = temp_TF
  final_effect_snp_frame[i,]$Motif_and_TF = temp_TF_motif
  final_effect_snp_frame[i,]$FDR_corr_P_value = temp_pval_corr
  final_effect_snp_frame[i,]$Raw_P_value = temp_pval
}

Ci_effect_SNPs = final_effect_snp_frame %>% 
  distinct() %>% 
  group_by(TF, SNP) %>%
  filter(FDR_corr_P_value == min(FDR_corr_P_value)) %>%
  ungroup()

write.table(Ci_effect_SNPs, "Ci_effect_SNPs.txt", row.names = FALSE, sep = "\t")

log_lik_ratio <- Ci_effect_SNPs %>%
  left_join(results_pval_sig, by = c("SNP" = "snpid", "Motif_and_TF" = "motif")) %>%
  dplyr::select(SNP, Motif_and_TF, log_lik_ratio) %>% distinct()

log_lik_ratio$TFSNP = paste(sapply(strsplit(log_lik_ratio$Motif_and_TF, "\\."), "[", 3), log_lik_ratio$SNP, sep = "-")

# taking the sum of LLR for the same snp-tf pairs with different motifs
log_lik_ratio <- log_lik_ratio %>%
  group_by(TFSNP) %>%
  mutate(total_log_lik_ratio = sum(log_lik_ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-log_lik_ratio, -Motif_and_TF) %>%
  distinct()

colnames(log_lik_ratio)[colnames(log_lik_ratio) == "total_log_lik_ratio"] <- "log_lik_ratio"
log_lik_ratio = tibble::column_to_rownames(log_lik_ratio, var = "TFSNP")

write.table(log_lik_ratio, "log_lik_ratio.txt", sep = "\t")
