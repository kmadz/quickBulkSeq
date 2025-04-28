#' quickHeatmap
#'
#' @param file A DESeqDataSet or path to a DESeq .rds object
#' @param title Character, title of plot
#' @param targets Character list, list of target genes to construct heatmap from, defaults to top 50 if empty
#' @param info Basis of differential expression analysis, defaults to "condition"
#' @param clustRows Logical, whether to cluster rows or not. default TRUE
#' @param clustCols Logical, whether to cluster columns or not, default TRUE
#' @param labelRows Logical, whether to label rows or not, default TRUE
#' @param n Integer, number of rows to display, default 50
#' @param num Character, the numerator of the comparison
#' @param denom Character, the denominator of the comparison
#' @param color A color palette, default colorRampPalette(list.reverse(brewer.pal(11, "PRGn")))(100)
#'
#' @return A heatmap
#' @export
#'
quickHeatmap <- function(file = NULL,
                          title = "",
                          targets = NULL,
                          info = 'condition',
                          clustRows = TRUE,
                          clustCols = TRUE,
                          labelRows = TRUE,
                          n = 50,
                          num = "dcas9_ANRIL",
                          denom = "dcas9_NONE",
                          color = colorRampPalette(list.reverse(brewer.pal(11, "PRGn")))(100)) {

  ### replace here the file handler for .rds vs dataframe



  # if (is.character(dds)) {
  #
  #   if (saveToDir) {
  #     output <- dirname(dds_path)
  #     title <- gsub(".*dds_", "", dds_path)
  #     title <- gsub(".rds", "", title)
  #     if (grepl("_", title)) {
  #       title <- gsub("_", " ", title)
  #     }
  #   }
  # } else if (!inherits(dds, "DESeqDataSet")) {
  #   stop("dds must be either a DESeqDataSet object or a valid file path to an RDS object.")
  # }

  dir <- dirname(file)
  dds <- readRDS(file)

  labs <- gsub(".*/results/", "", dir)
  if (grepl("_", labs)) {
    labs <- gsub("_", " ", labs)
  }

  # save the colData
  annot_info <- as.data.frame(colData(dds)[info])

  # save a Variance Stabilized Transformation of dds
  vsd <- vst(dds, blind = FALSE)

  res <- lfcShrink(
    dds,
    alpha = 0.05,
    contrast = c("condition", num, denom),
    parallel = FALSE,
    type = "ashr"
  )

  # Convert results to data frame and add gene_id
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    left_join(g2name, by = "gene_id") %>%
    arrange(pvalue)

  # Use normalized counts
  norm_df <- as.data.frame(counts(dds, normalized = TRUE)) %>%
    rownames_to_column("gene_id")

  # Merge results with assay data
  res <- left_join(res_df, norm_df, by = "gene_id")

  filtered_res <- na.omit(res)
  filtered_res <- filtered_res[filtered_res$padj < 0.05, ]

  if (length(filtered_res$padj) <= 1) {
    message("no")
    return()
  }
  # top_hits_fifty <- filtered_res[order(filtered_res$log2FoldChange), ][1:n, ]
  # targets <- rownames(top_hits_fifty)
  # targets <- rownames(filtered_res)

  targets <- filtered_res$gene_id
  numSigs <- length(targets)
  selected_genes <- g2name %>%
    filter(gene_id %in% targets) %>%
    arrange(factor(gene_name, levels = targets)) %>%
    distinct(gene_name, .keep_all = TRUE)


  # Check which gene_ids from selected_genes are found in the rownames of vsd
  matching_genes <- intersect(selected_genes$gene_id, rownames(assay(vsd)))

  # Subset the assay data using the correct rownames (matching_genes)
  selected_counts <- assay(vsd)[match(matching_genes, rownames(assay(vsd))), ]


  # define z score function
  zscore <- t(scale(t(selected_counts)))

  if (save) {
    png(
      paste(dir, "/", title, " ", labs, "Heatmap.png", sep = ""),
      width = 2400,
      height = 1800,
      res = 300
    )
  }

  pheatmap(
    zscore,
    color = color,
    cluster_rows = clustRows,
    show_rownames = labelRows,
    labels_row = selected_genes$gene_name,
    cluster_cols = clustCols,
    border_color = NA,
    annotation_col = annot_info,
    treeheight_row = 0,
    treeheight_col = 0,
    #annotation_colors = ann_colors,
    main = paste(labs, " ", title, ", ", numSigs, " significant genes", sep = "")  # Title of the heatmap
  )

  if (save) {
    dev.off()
  }

}

