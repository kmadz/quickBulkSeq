#' quickHeatmap
#'
#' @param file A DESeqDataSet or path to a DESeq .rds object
#' @param title Character, title of plot
#' @param targets Character list, list of target genes to construct heatmap from, defaults to top 50 if empty
#' @param info Basis of differential expression analysis, defaults to "condition"
#' @param clustRows Logical, whether to cluster rows or not. default TRUE
#' @param clustCols Logical, whether to cluster columns or not, default TRUE
#' @param labelRows Logical, whether to label rows or not, default TRUE
#' @param numRows Integer, number of rows to display, default 50
#' @param num Character, the numerator of the comparison
#' @param denom Character, the denominator of the comparison
#' @param color A color palette, default colorRampPalette(list.reverse(brewer.pal(11, "PRGn")))(100)
#' @param output Character, path to output directory, default current directory
#' @param padj Double, value for p-adjusted value threshold, default 0.05
#' @param save Logical, whether to save to output directory or not, default TRUE
#'
#' @importFrom pheatmap pheatmap
#' @importFrom DESeq2 DESeq
#' @importFrom SummarizedExperiment colData assay
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlist list.reverse
#'
#' @import dplyr
#' @return A heatmap
#' @export
#'
quickHeatmap <- function(file = NULL,
                         title = "",
                         output = ".",
                         targets = NULL,
                         info = 'condition',
                         clustRows = TRUE,
                         clustCols = TRUE,
                         labelRows = TRUE,
                         numRows = 50,
                         num = "",
                         denom = "",
                         padj = 0.05,
                         save = TRUE,
                         color = colorRampPalette(list.reverse(brewer.pal(11, "PRGn")))(100)) {
  if (num == "" | denom == "") {
    stop("Numerator and denominator for results not specified")
  }

  input <- ddsHandler(file, ".", title = title)
  annot_info <- as.data.frame(colData(input$dds)[info])


  # save a Variance Stabilized Transformation of dds
  vsd <- vst(input$dds, blind = FALSE)

  res <- lfcShrink(
    input$dds,
    alpha = 0.05,
    contrast = c("condition", num, denom),
    parallel = FALSE,
    type = "ashr"
  )

  # converts results to tibble and arrange by p value
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    left_join(g2name, by = "gene_id") %>%
    arrange(padj)

  # get normalized counts
  norm_df <- as.data.frame(counts(input$dds, normalized = TRUE)) %>%
    rownames_to_column("gene_id")

  # merge normalized counts to results
  res <- left_join(res_df, norm_df, by = "gene_id")

  filtered_res <- na.omit(res)
  filtered_res <- filtered_res[filtered_res$padj < padj, ]


  if (length(filtered_res$padj) <= 1) {
    message("no")
    return()
  }

  if (is.null(targets)) {
    targets <- filtered_res$gene_name[1:numRows]
  }

  numSigs <- nrow(filtered_res)


  selected_genes <- g2name %>%
    filter(gene_name %in% targets) %>%
    arrange(factor(gene_name, levels = targets)) %>%
    distinct(gene_name, .keep_all = TRUE)


  matching_genes <- intersect(selected_genes$gene_id, rownames(assay(vsd)))
  selected_counts <- assay(vsd)[match(matching_genes, rownames(assay(vsd))), ]


  # define z score function
  zscore <- t(scale(t(selected_counts)))

  png(
    file.path(input$output, paste(input$title, "Heatmap.png", sep = "_")),
    width = 2400,
    height = 1800,
    res = 300,
    units = "px"
  )

  final <- pheatmap(
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
    main = paste(input$title, "|", numSigs, "significant genes", sep = " ")  # Title of the heatmap
  )

  dev.off()

  return (final)
}

