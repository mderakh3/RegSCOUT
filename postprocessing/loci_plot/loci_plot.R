setwd("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/loci_plot/")

library(ggrepel)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(Seurat, warn.conflicts = FALSE)
library(Signac, warn.conflicts = FALSE)
library(readxl)
library(ape)
library(ggpubr)
library(fastmatch)
library(png)
library(grid)
library(patchwork)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(tibble)
library(Matrix)

# source functions
source("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/loci_plot/PlotSourceFilt.R")

cell_peaks <- read_excel("journal.pgen.1010759.s019.xlsx", sheet = 1)
full_table <- read.csv("complete_table.csv")
prio_table <- read.table("prioritized_table.txt", header = TRUE, sep = "\t")
loci_head <- read_excel("fine_loci_head.xlsx", sheet = 1)
eff_snp <- read.table("Ci_effect_SNPs.txt", header = TRUE, sep = "\t")
eff_snp <- eff_snp %>% dplyr::select(SNP, CHR, Pos) %>% distinct()

# generate a peaks column for next for loop
cell_peaks$peak <- paste(cell_peaks$chr, cell_peaks$start, cell_peaks$end, sep = "-")

# make a list of peaks for each cell type
cell_peaks_sep <- cell_peaks %>%
  separate_rows(`cell sub-types`, sep = ",")

peaks_grouped_by_cell <- cell_peaks_sep %>%
  group_by(`cell sub-types`) %>%
  summarize(peaks = list(unique(peak)), .groups = 'drop')

peaks_grouped_list <- split(peaks_grouped_by_cell$peaks, peaks_grouped_by_cell$`cell sub-types`)

# creating the Seurat object for snATAC
pbmc_matrix <- cell_peaks_sep %>%
  dplyr::select(-chr, -start, -end) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = `cell sub-types`, values_from = value, values_fill = list(value = 0)) %>%
  column_to_rownames(var = "peak") %>%
  as.matrix() %>%
  as("dgCMatrix")

chrom_assay <- CreateChromatinAssay(
  counts = pbmc_matrix,
  sep = c("-", "-"),
  genome = 'hg19',  # Change this to your genome build
  fragments = NULL, # Add fragment file if available
  min.cells = 1,
  min.features = 1
)

pbmc <- CreateSeuratObject(counts = pbmc_matrix, project = "96k")
pbmc[["peaks"]] <- chrom_assay
DefaultAssay(pbmc) <- "peaks"

# add effect snp position to full_table
full_table$effect_snp_pos = NA
for (i in 1:nrow(eff_snp)){
  snp = eff_snp$SNP[i]
  chr = eff_snp$CHR[i]
  pos = eff_snp$Pos[i]
  full_table$Pos[full_table$effect_snp == snp] = pos
  full_table$CHR[full_table$effect_snp == snp] = chr
}

loc_id = 129
loc_table = full_table[full_table$locus == loc_id,]
prio_loc_table = prio_table[prio_table$locus == loc_id,]
locus_region = strsplit(loc_table$loci_region[1], '-')[[1]]
loc_chr = locus_region[1]
loc_start = as.integer(locus_region[2])
loc_end = as.integer(locus_region[3])

peak_list = paste(loc_table[loc_table$rmp != 0,]$rmp, 
                  collapse = ",")
peak_list = str_split(peak_list, ",")[[1]]
peak_list = unique(peak_list)

peak_granges = StringToGRanges(peak_list)

# defining a zoomed in region and a main region for locus plot
zoom_reg = paste(seqnames(peak_granges)[1], (min(start(peak_granges)) - 50000),
                 (max(end(peak_granges)) + 50000), sep = "-")
middle_of_zoom <- round(min(start(peak_granges)) + (max(end(peak_granges)) - min(start(peak_granges)))/2)

# defining the main region to capture all relevant info including genes prioritized
# and to be centred around zoom_reg
interacting_promoters <- unique(loc_table$cicero_interacting_prom)
interacting_promoters <- Filter(function(x) x != "-", interacting_promoters)
temp_regions_df <- data.frame(peaks = c(unique(prio_loc_table$loci_region),unique(loc_table$rmp),interacting_promoters))
temp_regions_df <- temp_regions_df %>% distinct() %>%
  separate(peaks, into = c("chr", "start", "end"), sep = "-")

start_dist_from_mid <- middle_of_zoom - min(as.integer(temp_regions_df$start))
end_dist_from_mid <- max(as.integer(temp_regions_df$end)) - middle_of_zoom

if (start_dist_from_mid > end_dist_from_mid) {
  main_reg = paste(temp_regions_df$chr[1], min(as.integer(temp_regions_df$start)) - 50000, middle_of_zoom + start_dist_from_mid + 50000, sep = '-')
} else {
  main_reg = paste(temp_regions_df$chr[1], middle_of_zoom - end_dist_from_mid - 50000, max(as.integer(temp_regions_df$end)) + 50000, sep = '-')
}

cell_types = unique(full_table$cell_type)

# creating a list of chrs, start pos, and end pos for creation of peak plot
list_of_chrs = list()
list_of_starts = list()
list_of_ends = list()

# the desired region in which peaks must fall in to be included in peak plot
desired_chr = as.character(str_split(zoom_reg, '-')[[1]][1])
desired_start = as.integer(str_split(zoom_reg, '-')[[1]][2])
desired_end = as.integer(str_split(zoom_reg, '-')[[1]][3])

# loop to store all open peaks that fall into the zoom-reg
for (cell in names(peaks_grouped_list)){
  print(cell)
  cell_regions = peaks_grouped_list[[cell]][[1]]
  temp_starts_list = list()
  temp_ends_list = list()
  temp_chr_list = list()
  for (peak in cell_regions) {
    test = c(str_split(peak, '-')[[1]][1], str_split(peak, '-')[[1]][2], str_split(peak, '-')[[1]][3])
    if (test[1] == desired_chr & as.integer(test[2]) > desired_start & as.integer(test[3]) < desired_end) {
      print('success')
      temp_chr_list[length(temp_chr_list)+1] = as.character(test[1])
      temp_starts_list[length(temp_starts_list)+1] = as.integer(test[2])
      temp_ends_list[length(temp_ends_list)+1] = as.integer(test[3])
    }
  }
  list_of_chrs[length(list_of_chrs)+1] = list(temp_chr_list)
  list_of_starts[length(list_of_starts)+1] = list(temp_starts_list)
  list_of_ends[length(list_of_ends)+1] = list(temp_ends_list)
}

yaxisrange <- seq(1, length(cell_types))
xaxisrange <- seq(desired_start, desired_end, length.out = length(yaxisrange))
yaxislabels <- cell_types
list_of_diff_colours <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                          "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                          "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                          "#8A7C64", "#599861")

df <- data.frame(xaxisrange, yaxisrange)

# to ensure a different colour is picked for each cell type
segments_df <- data.frame()
number_of_colours = 20
for (i in 1:length(yaxislabels)) {
  if (length(list_of_starts[[i]]) == 0) {
    next
  } else {
    print(i)
    number_for_colour <- sample(1:number_of_colours, 1)
    peak_colour <- list_of_diff_colours[number_for_colour]
    list_of_diff_colours <- list_of_diff_colours[list_of_diff_colours != peak_colour]
    start_peaks <- as.integer(list_of_starts[[i]])
    end_peaks <- as.integer(list_of_ends[[i]])
    
    if (length(start_peaks) == 0) {
      print('no peaks to plot for this cell type')
    } else {
      temp_df <- data.frame(
        start = start_peaks,
        end = end_peaks,
        y = i,
        color = peak_colour
      )
      segments_df <- rbind(segments_df, temp_df)
    }
    number_of_colours = number_of_colours - 1
  }
}

# making the peak plot (maybe make it so there is a legend that says blue means locus)
peak_plot_zoom <- ggplot(df, aes(xaxisrange, yaxisrange)) + geom_blank() +
  scale_x_continuous(limits=c(desired_start, desired_end)) +
  scale_y_discrete(limits = factor(seq(1, length(yaxislabels), by = 1)), labels = yaxislabels) + 
  xlab(paste0(as.character(desired_chr), " position (bp)")) +
  ylab("Cell Type") +
  theme_linedraw() +
  geom_segment(data = segments_df, aes(x = start, y = y, xend = end, yend = y, color = color), linewidth = 2) +
  geom_vline(xintercept = unique(loc_table$Pos), linetype = "dashed", color = "red", size = 0.5) +
  theme(legend.position = 'none')

if (loc_start <= min(start(peak_granges)) - 50000) {
  peak_plot_zoom = peak_plot_zoom + annotate("rect", xmin = -Inf, xmax = loc_end, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightblue")
} else if (loc_end >= max(end(peak_granges)) + 50000) {
  peak_plot_zoom = peak_plot_zoom + annotate("rect", xmin = loc_start, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightblue")
} else {
  peak_plot_zoom = peak_plot_zoom + annotate("rect", xmin = loc_start, xmax = loc_end, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightblue")
}

# set positions of SNPs in next plot, modified x_pos to spread SNP names out more
loc_effect_snp <- loc_table[,c('effect_snp','CHR','Pos','tf')]
loc_effect_snp <- loc_effect_snp %>% distinct()
loc_effect_snp <- loc_effect_snp[order(loc_effect_snp$Pos), ]
snp_pos = start_pos_det(loc_effect_snp$Pos, rep(0, length(loc_effect_snp$Pos)))

snp_x_pos = snp_pos[[1]]
snp_y_pos = snp_pos[[2]]

# to space SNP RSID labels farther apart in locus plot
adjust_x_positions <- function(x_positions, spacing = 10000) {
  n <- length(x_positions)
  adjusted_positions <- x_positions
  if (n > 1) {
    for (i in 2:n) {
      if (adjusted_positions[i] - adjusted_positions[i - 1] < spacing) {
        adjusted_positions[i] <- adjusted_positions[i - 1] + spacing
        adjusted_positions[i-1] <- adjusted_positions[i - 1] - spacing
      }
    }
  }
  return(adjusted_positions)
}

snp_x_pos_mod = adjust_x_positions(loc_effect_snp$Pos)

snp_ranges = IRanges(start = snp_x_pos, end = (snp_x_pos+1))
snp_granges = GRanges(seqnames = loc_effect_snp$CHR,
                      ranges = snp_ranges)
snp_obj_count = matrix(0, nrow = length(snp_granges), ncol = ncol(pbmc))
row.names(snp_obj_count) = GRangesToString(snp_granges)
colnames(snp_obj_count) = colnames(pbmc)
snp_assay = CreateChromatinAssay(counts = snp_obj_count)
snp_obj = CreateSeuratObject(counts = snp_assay)

# making plot with SNP names
snp_name_plot = PeakPlot(
  object = snp_obj,
  region = zoom_reg
)+ annotate("text", x=snp_x_pos_mod, y=snp_y_pos,
            label= loc_effect_snp$effect_snp,
            angle = 90, size = 3.5) +
  annotate("segment", x = loc_effect_snp$Pos, y = -Inf, 
           xend = snp_x_pos_mod, yend = (snp_y_pos - 0.65),
           size = 0.7,
           color = "red")+
  ylim(c(-1,1)) #+ xlim(desired_start, desired_end)

snp_name_plot$labels$y = "RSID"
print(snp_name_plot)

snp_tf_frame = as.data.frame(matrix(0, nrow = nrow(loc_effect_snp),
                                    ncol = 4))
colnames(snp_tf_frame) = c("chr","start","end","TF_name")
snp_tf_frame$chr = loc_effect_snp$CHR
snp_tf_frame$start = loc_effect_snp$Pos - 10
snp_tf_frame$end = loc_effect_snp$Pos + 10
snp_tf_frame$TF_name = loc_effect_snp$tf
snp_tf_frame$TF_name = gsub(",","\n",snp_tf_frame$TF_name)

snp_tf_ranges = IRanges(start = snp_tf_frame$start,
                        end = snp_tf_frame$end)
snp_tf_granges = GRanges(seqnames = snp_tf_frame$chr,
                         ranges = snp_tf_ranges)

snp_tf_count = matrix(0, nrow = length(snp_tf_granges), ncol = ncol(pbmc))
colnames(snp_tf_count) = colnames(pbmc)
rownames(snp_tf_count) = GRangesToString(snp_tf_granges)
snp_tf_assay = CreateChromatinAssay(counts = snp_tf_count) 
snp_tf_obj = CreateSeuratObject(counts = snp_tf_assay)

snp_tf_plot <- PeakPlot(
  object = snp_obj,
  region = zoom_reg
) + annotate("text", x=snp_x_pos_mod, y=snp_y_pos, label= snp_tf_frame$TF_name, size = 2)+
  ylim(c(-1,1))

snp_tf_plot$labels$y = "Affected TFs"
print(snp_tf_plot)

cov_plot = CombineTracks(plotlist = list(snp_tf_plot, snp_name_plot, peak_plot_zoom),
                         heights = c(2,2,4), width = c(10,10,10))
cov_plot

reg = main_reg

# making GRange for reg
xnewtest = str_split(reg,"-", simplify = TRUE)
region_frame = data.frame(matrix(NA, nrow = 1, ncol=3))
region_frame$chr = xnewtest[1:nrow(xnewtest),1]
main_chr = xnewtest[1:nrow(xnewtest),1]
region_frame$start = xnewtest[1:nrow(xnewtest),2]
region_frame$end = xnewtest[1:nrow(xnewtest),3]
region_ranges = IRanges(start = as.integer(region_frame$start), 
                        end = as.integer(region_frame$end))
region_granges = GRanges(seqnames = region_frame$chr, ranges = region_ranges)

pbmc_regions = cell_peaks$peak

# making GRange for peaks in ATAC-seq data
xnewtest = str_split(pbmc_regions,"-", simplify = TRUE)
region_frame = data.frame(matrix(NA, nrow = length(pbmc_regions), ncol=3))
colnames(region_frame) = c("chr","start","end")
region_frame$chr = xnewtest[1:nrow(xnewtest),1]
main_chr = xnewtest[1:nrow(xnewtest),1]
region_frame$start = xnewtest[1:nrow(xnewtest),2]
region_frame$end = xnewtest[1:nrow(xnewtest),3]

peak_ranges = IRanges(start = as.integer(region_frame$start), 
                      end = as.integer(region_frame$end))
peak_granges = GRanges(seqnames = region_frame$chr, ranges = peak_ranges)


snp_ranges = IRanges(start = (eff_snp$Pos - 1), end = eff_snp$Pos)
snp_granges = GRanges(seqnames = eff_snp$CHR, ranges = snp_ranges)

snp_peak_overlap = findOverlaps(snp_granges, peak_granges)
snp_index = queryHits(snp_peak_overlap)
snp_granges = unique(snp_granges[snp_index])

tf_ranges = IRanges(start = (loc_effect_snp$Pos - 500), end = (loc_effect_snp$Pos + 500))
tf_granges = GRanges(seqnames = loc_effect_snp$CHR, ranges = tf_ranges)
tf_granges$TF_id = loc_effect_snp$tf

tf_peak_overlap = findOverlaps(tf_granges, peak_granges)
tf_index = queryHits(tf_peak_overlap)
tf_granges = unique(tf_granges[tf_index])

reg_snp_overlap = findOverlaps(region_granges, snp_granges)
snp_list = snp_granges[unique(subjectHits(reg_snp_overlap)),]

tf_snp_overlap = findOverlaps(region_granges, tf_granges)
tf_list = tf_granges[unique(subjectHits(tf_snp_overlap))]

tf_frame = as.data.frame(matrix(nrow = length(tf_list), ncol = 2))
colnames(tf_frame) = c("TF_name", "start")

tf_frame$TF_name = tf_list$TF_id
tf_frame$start = start(tf_list)

# for peak_plot
snp_y_start = c(0.2,0.2)

tf_frame = tf_frame[order(tf_frame$start),]
tf_number = nrow(tf_frame)
y_tf_step = seq(0.1, (tf_number*0.1), by = 0.1)
x_tf_start = seq((min(tf_frame$start) - 5e4),(max(tf_frame$start) + 5e4)
                 , by=((max(tf_frame$start)-min(tf_frame$start)+1e5)/(tf_number - 1)))

peak_plot = PeakPlot(
  object = pbmc,
  region = reg
) + ylim(-1,1)
peak_plot$labels$y = "Peaks"

main_reg_start = as.integer(str_split(main_reg, '-')[[1]][2])
main_reg_end = as.integer(str_split(main_reg, '-')[[1]][3])
peak_plot_starts = as.integer(c(desired_start, desired_end))
peak_plot_ends = as.integer(c(main_reg_start+10000, main_reg_end-10000))
peak_plot = peak_plot + annotate("segment", x = peak_plot_starts, y = snp_y_start,
                                 xend = peak_plot_ends, yend = c(+Inf,+Inf), 
                                 color = 'black', size = 1) +
  geom_vline(xintercept=loc_effect_snp$Pos, linetype="dashed", 
             color = "red", size=0.5)

peak_plot = peak_plot + annotate("rect", xmin = loc_start, xmax = loc_end, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightblue")
peak_plot

# defining prioritized genes
prio_genes <- unique(prio_loc_table$gene_score)
prio_genes <- gsub("\\(.*?\\)", "", prio_genes)
prio_genes <- gsub(" ", "", prio_genes)
prio_genes <- str_split(prio_genes, ",", simplify = TRUE)
prio_genes <- prio_genes[prio_genes != ""]
prio_genes <- unique(prio_genes)

# reading in gencode
genecode_dir = "/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/loci_plot/gencode.v19.annotation.gff3"
gene_code_data = read.gff(genecode_dir, na.strings = c(".", "?"), GFF3 = TRUE)
gene_code_data = subset(gene_code_data, type != "gene")

new_annotation_ranges = IRanges(start = gene_code_data$start,
                                end = gene_code_data$end)

new_annotation = GRanges(ranges = new_annotation_ranges, 
                         seqnames = gene_code_data$seqid, 
                         strand = gene_code_data$strand)

new_annotation$type = gene_code_data$type

gene_name_list = str_split(gene_code_data$attributes, "gene_name=", simplify = TRUE)
gene_name_list = str_split(gene_name_list[,2], ";", simplify = TRUE)[,1]
print(head(gene_name_list,10))

gene_type_list = str_split(gene_code_data$attributes, "gene_type=", simplify = TRUE)
gene_type_list = str_split(gene_type_list[,2], ";", simplify = TRUE)[,1]
print(head(gene_type_list,10))

gene_id_list = str_split(gene_code_data$attributes, "gene_id=", simplify = TRUE)
gene_id_list = str_split(gene_id_list[,2], ";", simplify = TRUE)[,1]
print(head(gene_id_list,10))

tx_id_list = str_split(gene_code_data$attributes, "transcript_id=", simplify = TRUE)
tx_id_list = str_split(tx_id_list[,2], ";", simplify = TRUE)[,1]
print(head(tx_id_list,10))

new_annotation$gene_name = gene_name_list
new_annotation$gene_biotype = gene_type_list
new_annotation$gene_id = gene_id_list
new_annotation$tx_id = tx_id_list

# make annotation frame object, can either use annotation from pbmc or the gencode you read in
annotation_frame = new_annotation

#To make annotation_frame_all protein coding specific
procode_index_2 = annotation_frame$gene_biotype == "protein_coding"
annotation_frame = annotation_frame[procode_index_2]

reg_granges_all = StringToGRanges(reg)
reg_annot = findOverlaps(reg_granges_all, annotation_frame)
found_gene = annotation_frame[unique(subjectHits(reg_annot))]
all_gene_number = length(unique(found_gene$gene_name))

# to show all genes instead of only significantly associated genes
Annotation(pbmc) <- annotation_frame

gene_plot <- AnnotationPlotRedLine(  
  object = pbmc,
  region = reg
)

Ypos = (gene_plot$plot_env$y_limit[2] - 
          gene_plot$plot_env$y_limit[1])/2 + gene_plot$plot_env$y_limit[1]

gene_plot <- gene_plot + annotate("rect", xmin = loc_start, xmax = loc_end, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "lightblue")

second_half_of_plot = CombineTracks(
  plotlist = list(peak_plot, gene_plot),
  height = c(1,4)
)

# using this function keeps the X-axis labels in the plot
main_plot = CombineTracksKeepXAxis(
  plotlist = list(cov_plot, second_half_of_plot),
  heights = c(3,3)
)

loc_fig_dir = paste0("/Users/m.hossein_drn/Documents/aim_1/Manuscript/tables_figures/loci_plot/",
                     loc_id,".png")

png(filename = loc_fig_dir, width = 3600, height = 4500, res = 300)

main_plot

dev.off()

