setwd("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/part_1")
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

fine_loci_head <- readxl::read_excel("fine_loci_head.xlsx")

UC_loci <- fine_loci_head[fine_loci_head$Trait %in% c("UC"), ]
CD_loci <- fine_loci_head[fine_loci_head$Trait %in% c("CD"), ]
IBD_loci <- fine_loci_head[fine_loci_head$Trait %in% c("IBD"), ]

ci_snps <- read.table("ci_snp_locus.txt", header = T)

# CI-SNPs distribution in traits plot
CD = length(ci_snps$id[ci_snps$trait == "CD"])
UC = length(ci_snps$id[ci_snps$trait == "UC"])
IBD = length(ci_snps$id[ci_snps$trait == "IBD"])

data <- data.frame(
  Trait = c("CD", "UC", "IBD"),
  Count = c(CD, UC, IBD))

p <- ggplot(data, aes(x = Trait, y = Count, fill = Trait)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("blue", "#00B050", "red")) +
  labs(y = "Number of CI-SNPs") +
  theme_minimal(base_size = 25) +
  theme(
    text = element_text(color = "black"),
    #plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 30, margin = margin(r = 20)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25, color = "black"),
    axis.text.y = element_text(size = 30, color = "black"),
    axis.title.x = element_blank(),
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    #plot.background = element_rect(color = "black", linewidth = 1.5)
  )

print(p)
ggsave("CI-SNPs_traits.png", plot = p, width = 10, height = 7, dpi = 300)

# CI-SNPs distribution in studies plot
GCST004131 <- length(ci_snps$id[ci_snps$author == "delange"])
GCST005837  <- length(ci_snps$id[ci_snps$author == "mark"])
Merged <- length(ci_snps$id)

data <- data.frame(
  Study = c("GCST004131", "GCST005837", "Merged"),
  Count = c(GCST004131, GCST005837, Merged))

p <- ggplot(data, aes(x = Study, y = Count, fill = Study)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("blue", "#00B050", "red")) +
  labs(y = "Number of CI-SNPs") +
  theme_minimal(base_size = 25) +
  theme(
    text = element_text(color = "black"),
    #plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 25, margin = margin(r = 20)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black"),
    axis.title.x = element_blank(),
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    #plot.background = element_rect(color = "black", linewidth = 1.5)
  )

print(p)
ggsave("CI-SNPs_studies.png", plot = p, width = 10, height = 7, dpi = 300)

# PPA-Locus for CI-SNPs scatter plot
colnames(ci_snps)[colnames(ci_snps) == "trait"] <- "Trait"

p <- ggplot(ci_snps, aes(x = locus, y = PPA, color = Trait)) +
  geom_point(size = 4) +  # Increase point size
  scale_color_manual(values = c("blue", "#00B050", "red")) +
  #labs(title = "PPA-Locus for CI-SNPs") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black", size = 1) +  # Add threshold line
  annotate("text", x = max(ci_snps$locus) - 10, y = 0.9, label = "PPA = 0.9", color = "black", size = 6, hjust = 1, vjust = -1) +  # Add label for threshold line
  theme_minimal(base_size = 30) +
  theme(
    text = element_text(color = "black"),
    #plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 30, margin = margin(r = 20)),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25, color = "black"),  # Rotate x-axis labels and adjust size
    axis.text.y = element_text(size = 25, color = "black"),
    axis.title.x = element_text(size = 25, margin = margin(t = 20)),
    plot.background = element_rect(fill = "white", color = NA)  # Set background to white
  ) +
  scale_x_continuous(breaks = seq(0, max(ci_snps$locus), by = 10))  # Set x-axis labels at intervals of 10

print(p)
ggsave("PPA_Locus_CI_SNPs.png", plot = p, width = 24, height = 10)  # Adjust width and height as needed

# A scatter plot of PPA-chr for CI-SNPs for IBD, CD, and UC
ci_snps$chr <- factor(ci_snps$chr, levels = paste0("chr", 1:22))

p <- ggplot(ci_snps, aes(x = chr, y = PPA, color = Trait)) +
  geom_point(size = 4, position = position_jitter(width = 0.3)) +  # Increase point size and add horizontal jitter
  scale_color_manual(values = c("blue", "#00B050", "red")) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black", size = 1) +  # Add threshold line
  annotate("text", x = 22.5, y = 0.9, label = "PPA = 0.9", color = "black", size = 6, hjust = -2, vjust = 2) +  # Add label for threshold line
  #labs(title = "PPA by Chromosome for CI-SNPs", x = "Chromosome", y = "PPA") +
  theme_minimal(base_size = 30) +
  theme(
    text = element_text(color = "black"),
    axis.title.y = element_text(size = 30, margin = margin(r = 20)),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 25, color = "black"),  # Rotate x-axis labels and adjust size
    axis.text.y = element_text(size = 25, color = "black"),
    axis.title.x = element_text(size = 25, margin = margin(t = 20)),
    plot.background = element_rect(fill = "white", color = NA)  # Set background to white
  )
print(p)

ggsave("PPA_by_Chromosome_CI_SNPs.png", plot = p, width = 24, height = 10)  # Adjust width and height as needed

ci_snps_table <- ci_snps %>%
  dplyr::select(id, PPA, chr, locus, Trait) %>%
  arrange(chr)

ci_snps_table <- ci_snps_table[ci_snps_table$PPA >= 0.9, ]
write.csv(ci_snps_table, "CI_SNPs_high_ppa.csv", row.names = FALSE)

# distribution plot for tfs
log_lik_ratio <- read.table("log_lik_ratio.txt", header = TRUE, sep = "\t")
log_lik_ratio <- tibble::rownames_to_column(log_lik_ratio, "TFSNP")

log_lik_ratio <- log_lik_ratio %>%
  separate(TFSNP, into = c("TF", "SNP"), sep = "-(?!.*-)", extra = "merge") %>% 
  dplyr::select(TF, SNP, log_lik_ratio)

effect_snp_count <- log_lik_ratio %>%
  group_by(TF) %>%
  summarize(count = n())

total_TFs_above_threshold <- effect_snp_count %>%
  filter(count >= 10) %>%
  nrow()

ggplot(effect_snp_count, aes(x = count)) +
  geom_line(stat = "count", color = "black", size = 1) +
  geom_point(stat = "count", color = "red", size = 3) +
  geom_area(stat = "count", aes(y = after_stat(count)), fill = "blue", alpha = 0.5, data = filter(effect_snp_count, count >= 10)) +
  labs(x = "Number of Effect-SNPs per TF", y = "Number of TFs") +
  theme_minimal(base_size = 35) + # Increase the base font size
  theme(
    axis.text = element_text(size = 35, color = "black"), # Increase font size and set color for axis text
    axis.title = element_text(size = 40, color = "black"), # Increase font size and set color for axis titles
    plot.title = element_text(size = 35, color = "black"), # Increase font size and set color for plot title
    legend.text = element_text(size = 35, color = "black"), # Increase font size and set color for legend text
    legend.title = element_text(size = 35, color = "black"), # Increase font size and set color for legend title
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white") # Set plot background to white
  ) +
  annotate("text", x = 10, y = 75, label = paste("Total TFs with >= 10 Effect-SNPs:", total_TFs_above_threshold), vjust = -1, hjust = 0.1, size = 12, color = "black") +
  annotate("segment", x = 10, xend = 10, y = 0, yend = 75, arrow = arrow(length = unit(0.2, "inches")), color = "black") +
  geom_hline(yintercept = 0, color = "black", size = 0.5) + # Add a solid black line on the x-axis
  ylim(0, 80)

ggsave("tfs_distribution.png", width = 20, height = 18)

# distribution plot for effect-snps
tf_count <- log_lik_ratio %>%
  group_by(SNP) %>%
  summarize(count = n())

total_SNPs_above_threshold <- tf_count %>%
  filter(count >= 10) %>%
  nrow()

ggplot(tf_count, aes(x = count)) +
  geom_line(stat = "count", color = "black", size = 1) +
  geom_point(stat = "count", color = "red", size = 3) +
  geom_area(stat = "count", aes(y = after_stat(count)), fill = "blue", alpha = 0.5, data = filter(tf_count, count >= 10)) +
  labs(x = "Number of TFs per Effect-SNP", y = "Number of Effect-SNPs") +
  theme_minimal(base_size = 35) + # Increase the base font size
  theme(
    axis.text = element_text(size = 35, color = "black"), # Increase font size and set color for axis text
    axis.title = element_text(size = 40, color = "black"), # Increase font size and set color for axis titles
    plot.title = element_text(size = 35, color = "black"), # Increase font size and set color for plot title
    legend.text = element_text(size = 35, color = "black"), # Increase font size and set color for legend text
    legend.title = element_text(size = 35, color = "black"), # Increase font size and set color for legend title
    panel.background = element_rect(fill = "white", color = "white"), # Set background to white
    plot.background = element_rect(fill = "white", color = "white") # Set plot background to white
  ) +
  annotate("text", x = 10, y = 1000, label = paste("Total Effect-SNPs with >= 10 TFs:", total_SNPs_above_threshold), vjust = 2, hjust = 0.1, size = 12, color = "black") +
  annotate("segment", x = 10, xend = 10, y = 0, yend = 800, arrow = arrow(length = unit(0.2, "inches")), color = "black") +
  geom_hline(yintercept = 0, color = "black", size = 0.5) # Add a solid black line on the x-axis
  
ggsave("effect_snps_distribution.png", width = 16, height = 12)

# effect-snps heatmap
library(tidyverse)

log_lik_ratio_matrix <- log_lik_ratio %>%
  spread(TF, log_lik_ratio) %>%
  column_to_rownames(var = "SNP") %>%
  as.matrix()

log_lik_ratio_matrix <- log_lik_ratio_matrix[rownames(log_lik_ratio_matrix) %in% tf_count$SNP[tf_count$count >= 10], 
                                             colnames(log_lik_ratio_matrix) %in% effect_snp_count$TF[effect_snp_count$count >= 10]]

log_lik_ratio_matrix <- log_lik_ratio_matrix %>% replace(is.na(.), 0) %>%
  as.data.frame() %>%
  select_if(function(x) any(x != 0)) %>%
  as.matrix()

library(ComplexHeatmap)
library(circlize)

log_lik_ratio_heatmap <- Heatmap(
  log_lik_ratio_matrix, 
  name = "LLR",
  col = circlize::colorRamp2(c(-20, 0, 20), c("blue", "white", "red")),
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(2, "cm"),
  border = TRUE, # Add borders around the entire heatmap
  rect_gp = gpar(col = "black", lwd = 1), # Add borders around the tiles
  row_title = "SNPs",
  column_title = "Transcription Factors",
  row_title_gp = gpar(fontsize = 15), # Increase font size of the y-axis title
  column_title_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 5)
)

png("log_lik_ratio_heatmap.png", width = 3600, height = 3600, res = 300)
draw(log_lik_ratio_heatmap)
dev.off()

# number of SNPs per locus plot
ci_snp_locus <- read.table("ci_snp_locus.txt", header = TRUE, sep = "\t")

ci_locus_count <- ci_snp_locus %>%
  group_by(locus) %>%
  summarize(count = n())

Ci_effect_SNPs <- read.table("Ci_effect_SNPs.txt", header = TRUE, sep = "\t")

Ci_effect_SNPs <- Ci_effect_SNPs %>%
  left_join(ci_snp_locus, by = c("SNP" = "id")) %>%
  select(SNP, locus) %>% distinct()

effect_locus_count <- Ci_effect_SNPs %>%
  group_by(locus) %>%
  summarize(count = n())

library(readxl)
risk_regions_ratio <- read_excel("risk_regions_ratio.xlsx")

risk_regions_ratio <- risk_regions_ratio %>%
  separate(TFSNP, into = c("TF", "SNP"), 
           sep = "-(?!.*-)", extra = "merge")

risk_regions_ratio <- risk_regions_ratio %>%
  left_join(ci_snp_locus, by = c("SNP" = "id")) %>%
  select(SNP, locus) %>% distinct()

po_locus_count <- risk_regions_ratio %>%
  group_by(locus) %>%
  summarize(count = n())

locus_count <- data.frame(locus = 1:160)

locus_count <- locus_count %>%
  left_join(ci_locus_count, by = c("locus" = "locus")) %>%
  left_join(effect_locus_count, by = c("locus" = "locus")) %>%
  left_join(po_locus_count, by = c("locus" = "locus")) %>%
  rename(ci = count.x, effect = count.y, po = count)

locus_count[is.na(locus_count)] <- 0

locus_count_long <- locus_count %>%
  gather(key = "SNP_type", value = "count", ci, effect, po)

locus_count_long$SNP_type <- gsub("ci", "CI-SNPs", locus_count_long$SNP_type)
locus_count_long$SNP_type <- gsub("effect", "Effect-SNPs", locus_count_long$SNP_type)
locus_count_long$SNP_type <- gsub("po", "PO-SNPs", locus_count_long$SNP_type)

ggplot(locus_count_long, aes(x = factor(locus), y = count, fill = SNP_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("CI-SNPs" = "#00B050", "Effect-SNPs" = "red", "PO-SNPs" = "blue")) +
  labs(x = "Locus", y = "Number of SNPs", fill = "SNP Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, color = "black", angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels
    axis.text.y = element_text(size = 25, color = "black"),
    axis.title = element_text(size = 30, color = "black"),
    legend.text = element_text(size = 25, color = "black"), 
    legend.title = element_text(size = 25, color = "black"),
    legend.position = c(0.1, 0.9), # Position the legend inside the plot at the top right
    panel.background = element_rect(fill = "white", color = "white"), 
    plot.background = element_rect(fill = "white", color = "white")
  )

ggsave("SNPs_per_locus.png", width = 45, height = 18)

# PO-SNPs table 
risk_regions_ratio <- read_excel("risk_regions_ratio.xlsx")

risk_regions_ratio <- risk_regions_ratio %>%
  separate(TFSNP, into = c("TF", "SNP"), sep = "-(?!.*-)", extra = "merge") %>%
  dplyr::select(SNP, region, cell) %>%
  distinct() 

write.table(risk_regions_ratio, "po_snps_tble.txt", sep = "\t", row.names = FALSE)

