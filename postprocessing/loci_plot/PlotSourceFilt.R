AnnotationPlot <- function(
    object,
    region,
    assay = NULL,
    mode = "gene",
    sep = c("-", "-"),
    extend.upstream = 0,
    extend.downstream = 0
) {
  if(mode == "gene") {
    collapse_transcript <- TRUE
    label <- "gene_name"
  } else if (mode == "transcript") {
    collapse_transcript <- FALSE
    label <- "tx_id"
  } else {
    stop("Unknown mode requested, choose either 'gene' or 'transcript'")
  }
  annotation <- Annotation(object = object)
  if (is.null(x = annotation)) {
    return(NULL)
  }
  region <- FindRegion(
    object = object,
    region = region,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  if (mode == "gene") {
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
      fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  } else {
    tx.keep <- unique(x = annotation.subset$tx_id)
    annotation.subset <- annotation[
      fmatch(x = annotation$tx_id, table = tx.keep, nomatch = 0L) > 0L
    ]
  }
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      collapse_transcript = collapse_transcript
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 3
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    
    if (n_stack == 1){
      annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.1
    }else if(n_stack == 2){
      annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.2
    }else{ 
      annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.4
    }
    
    p <- p + geom_text(
      data = annotation_df_list$labels,
      mapping = aes_string(x = "position", y = "dodge", label = label),
      size = 2.5
    )
    y_limit <- c(0.9, n_stack + 0.4)
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(start.pos, end.pos) +
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}

reformat_annotations <- function(
    annotation,
    start.pos,
    end.pos,
    collapse_transcript = TRUE
) {
  total.width <- end.pos - start.pos
  tick.freq <- total.width / 50
  annotation <- annotation[annotation$type == "exon"]
  exons <- as.data.frame(x = annotation)
  if (collapse_transcript) {
    annotation <- split(
      x = annotation,
      f = annotation$gene_name
    )
  } else {
    annotation <- split(
      x = annotation,
      f = annotation$tx_id
    )
  }
  annotation <- lapply(X = annotation, FUN = as.data.frame)
  
  # add gene total start / end
  gene_bodies <- list()
  for (i in seq_along(annotation)) {
    df <- data.frame(
      seqnames = annotation[[i]]$seqnames[[1]],
      start = min(annotation[[i]]$start),
      end = max(annotation[[i]]$end),
      strand = annotation[[i]]$strand[[1]],
      tx_id = annotation[[i]]$tx_id[[1]],
      gene_name = annotation[[i]]$gene_name[[1]],
      gene_biotype = annotation[[i]]$gene_biotype[[1]],
      type = "body"
    )
    # trim any that extend beyond region
    df$start <- ifelse(
      test = df$start < start.pos,
      yes = start.pos,
      no = df$start
    )
    df$end <- ifelse(
      test = df$end > end.pos,
      yes = end.pos,
      no = df$end
    )
    breaks <- split_body(df = df, width = tick.freq)
    df <- rbind(df, breaks)
    gene_bodies[[i]] <- df
  }
  gene_bodies <- do.call(what = rbind, args = gene_bodies)
  
  # record if genes overlap
  overlap_idx <- record_overlapping(
    annotation = gene_bodies,
    min.gapwidth = 1000,
    collapse_transcript = collapse_transcript
  )
  # overlap_idx <- overlap_idx
  if (collapse_transcript) {
    gene_bodies$dodge <- overlap_idx[gene_bodies$gene_name]
    exons$dodge <- overlap_idx[exons$gene_name]
  } else {
    gene_bodies$dodge <- overlap_idx[gene_bodies$tx_id]
    exons$dodge <- overlap_idx[exons$tx_id]
  }
  
  label_df <- gene_bodies[gene_bodies$type == "body", ]
  label_df$width <- label_df$end - label_df$start
  label_df$position <- label_df$start + (label_df$width / 2)
  
  onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
  onminus <- gene_bodies[gene_bodies$strand == "-", ]
  
  return(
    list(
      "labels" = label_df,
      "exons" = exons,
      "plus" = onplus,
      "minus" = onminus
    )
  )
}

record_overlapping <- function(
    annotation,
    min.gapwidth = 1000,
    collapse_transcript = TRUE
) {
  # convert back to granges
  annotation$strand <- "*"
  gr <- makeGRangesFromDataFrame(
    df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
  )
  # work out which ranges overlap
  collapsed <- reduce(
    x = gr, with.revmap = TRUE, min.gapwidth = 100000
  )$revmap
  idx <- seq_along(gr)
  for (i in seq_along(collapsed)) {
    mrg <- collapsed[[i]]
    for (j in seq_along(mrg)) {
      idx[[mrg[[j]]]] <- j
    }
  }
  if (collapse_transcript) {
    names(x = idx) <- gr$gene_name
  } else {
    names(x = idx) <- gr$tx_id
  }
  return(idx)
}

FindRegion <- function(
    object,
    region,
    sep = c("-", "-"),
    assay = NULL,
    extend.upstream = 0,
    extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}

split_body <- function(df, width = 1000) {
  wd <- df$end - df$start
  nbreak <- wd / width
  if (nbreak > 1) {
    steps <- 0:(nbreak)
    starts <- (width * steps) + df$start
    starts[starts > df$end] <- NULL
  } else {
    starts <- df$end
  }
  breaks <- data.frame(
    seqnames = df$seqnames[[1]],
    start = starts,
    end = starts + 1,
    strand = df$strand[[1]],
    tx_id = df$tx_id[[1]],
    gene_name = df$gene_name[[1]],
    gene_biotype = df$gene_biotype[[1]],
    type = "arrow"
  )
  return(breaks)
}

start_pos_det = function(x_data, y_data){
  x_tf_start = x_data
  y_tf_step = y_data
  
  x_tf_start_new = x_data
  y_tf_step_new = y_data
  
  jump_stat = FALSE
  
  for (i in c(1:(length(x_tf_start)-0))){
    if(jump_stat){
      jump_stat = FALSE
      next
    }
    temp_tf_start = x_tf_start[i]
    temp_tf_y = 0
    if (i == 1){
      if (abs(x_tf_start[i+1]-x_tf_start[i])<5500){
        temp_tf_start = x_tf_start[i] - 10000
        temp_tf_y = 0
        jump_stat = TRUE
      }
    }else if(i == length(x_tf_start)){
      if(abs(x_tf_start[i] - x_tf_start[i-1])<5500){
        temp_tf_start = x_tf_start[i] + 10000
        temp_tf_y = 0
        jump_stat = TRUE
      }
    }else{
      prev_stat = abs(x_tf_start[i] - x_tf_start[i-1])<5500
      next_stat = abs(x_tf_start[i+1]-x_tf_start[i])<5500
      if((prev_stat == TRUE) & (next_stat == TRUE)){
        if(prev_stat < next_stat){
          temp_tf_start = x_tf_start[i] + 5000
          temp_tf_y = 0.5
        }else{
          temp_tf_start = x_tf_start[i] - 5000
          temp_tf_y = 0.5
        }
        jump_stat = TRUE
      }else if(prev_stat == TRUE){
        temp_tf_start = x_tf_start[i] + 10000
        temp_tf_y = 0
        jump_stat = TRUE
      }else if(next_stat == TRUE){
        temp_tf_start = x_tf_start[i] - 10000
        temp_tf_y = 0
        jump_stat = TRUE
      }
    }
    x_tf_start[i] = temp_tf_start
    y_tf_step[i] = temp_tf_y 
    
  }
  
  return (list(x_tf_start, y_tf_step))
  
}


## Modifying a function from the signac package, credit to https://github.com/stuart-lab/signac
CombineTracksKeepXAxis <- function(
    plotlist,
    expression.plot = NULL,
    heights = NULL,
    widths = NULL
) {
  # remove any that are NULL
  nullplots <- sapply(X = plotlist, FUN = is.null)
  plotlist <- plotlist[!nullplots]
  heights <- heights[!nullplots]
  
  if (length(x = plotlist) == 1) {
    return(plotlist[[1]])
  }
  
  # remove x-axis from all but last plot
  # for (i in 1:(length(x = plotlist) - 1)) {
  #   plotlist[[i]] <- plotlist[[i]] + theme(
  #     axis.title.x = element_blank(),
  #     axis.text.x = element_blank(),
  #     axis.line.x.bottom = element_blank(),
  #     axis.ticks.x.bottom = element_blank()
  #   )
  # }
  
  # combine plots
  if (is.null(x = heights)) {
    # set height of first element to 10x more than other elements
    n.plots <- length(x = plotlist)
    heights <- c(8, rep(1, n.plots - 1))
  } else {
    if (length(x = heights) != length(x = plotlist)) {
      stop("Relative height must be supplied for each plot")
    }
  }
  if (!is.null(x = expression.plot)) {
    # align expression plot with the first element in plot list
    p <- (plotlist[[1]] + expression.plot) +
      plot_layout(widths = widths)
    
    n <- length(x = plotlist)
    heights.2 <- heights[2:n]
    p2 <- wrap_plots(plotlist[2:n], ncol = 1, heights = heights.2)
    
    p <- p + p2 + guide_area() + plot_layout(
      ncol = 2, heights = c(heights[[1]], sum(heights.2)),
      guides = "collect")
  } else {
    p <- wrap_plots(plotlist, ncol = 1, heights = heights)
  }
  return(p)
}

#making sure y-axis labels stay at y-axis even when plots are wrapped tgt
CombineTracksDefineMargins <- function(
    plotlist,
    expression.plot = NULL,
    heights = NULL,
    widths = NULL
) {
  # remove any that are NULL
  nullplots <- sapply(X = plotlist, FUN = is.null)
  plotlist <- plotlist[!nullplots]
  heights <- heights[!nullplots]
  
  if (length(x = plotlist) == 1) {
    return(plotlist[[1]])
  }
  
  # remove x-axis from all but last plot
  # for (i in 1:(length(x = plotlist) - 1)) {
  #   plotlist[[i]] <- plotlist[[i]] + theme(
  #     axis.title.x = element_blank(),
  #     axis.text.x = element_blank(),
  #     axis.line.x.bottom = element_blank(),
  #     axis.ticks.x.bottom = element_blank()
  #   )
  # }
  
  # combine plots
  if (is.null(x = heights)) {
    # set height of first element to 10x more than other elements
    n.plots <- length(x = plotlist)
    heights <- c(8, rep(1, n.plots - 1))
  } else {
    if (length(x = heights) != length(x = plotlist)) {
      stop("Relative height must be supplied for each plot")
    }
  }
  if (!is.null(x = expression.plot)) {
    # align expression plot with the first element in plot list
    p <- (plotlist[[1]] + expression.plot) +
      plot_layout(widths = widths)
    
    n <- length(x = plotlist)
    heights.2 <- heights[2:n]
    p2 <- wrap_plots(plotlist[2:n], ncol = 1, heights = heights.2)
    
    p <- p + p2 + guide_area() + plot_layout(
      ncol = 2, heights = c(heights[[1]], sum(heights.2)),
      guides = "collect")
  } else {
    p <- wrap_plots(plotlist, ncol = 1, heights = heights) & theme(plot.margin = margin(t = 0, r = 3, b = 0, l = 5))
  }
  return(p)
}
#& theme(plot.margin = margin(t = 0, r = 10, b = 0, l = 5))

##Also taken from Stuart lab, annotation plot but with a red line drawn before text
#Or else the red line would potentially cover gene names and gene diagrams
AnnotationPlotRedLine <- function(
    object,
    region,
    assay = NULL,
    mode = "gene",
    sep = c("-", "-"),
    extend.upstream = 0,
    extend.downstream = 0
) {
  if(mode == "gene") {
    collapse_transcript <- TRUE
    label <- "gene_name"
  } else if (mode == "transcript") {
    collapse_transcript <- FALSE
    label <- "tx_id"
  } else {
    stop("Unknown mode requested, choose either 'gene' or 'transcript'")
  }
  annotation <- Annotation(object = object)
  if (is.null(x = annotation)) {
    return(NULL)
  }
  region <- FindRegion(
    object = object,
    region = region,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  if (mode == "gene") {
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
      fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  } else {
    tx.keep <- unique(x = annotation.subset$tx_id)
    annotation.subset <- annotation[
      fmatch(x = annotation$tx_id, table = tx.keep, nomatch = 0L) > 0L
    ]
  }
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      collapse_transcript = collapse_transcript
    )
    p <- ggplot() +
      # My addition: dashed red line for where the effect SNPs are
      geom_vline(xintercept=loc_effect_snp$Pos, linetype="dashed", 
                 color = "red", size=0.5) +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 3
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.35  # was 0.2
    # print('here')
    # print(annotation_df_list$labels$gene_name %in% prio_genes)
    # print('here2')
    p <- p + geom_text(
      data = annotation_df_list$labels[!(annotation_df_list$labels$gene_name %in% prio_genes),],
      mapping = aes_string(x = "position", y = "dodge", label = label),
      size = 2.7, # used to be 2.5 (I changed it),
      color = 'black'
    ) + geom_text( # this geom_text is my addition as well
      data = annotation_df_list$labels[annotation_df_list$labels$gene_name %in% prio_genes,],
      mapping = aes_string(x = "position", y = "dodge", label = label),
      size = 2.7,
      color = 'red'
    )
    y_limit <- c(0.9, n_stack + 0.4)
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(start.pos, end.pos) +
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}
