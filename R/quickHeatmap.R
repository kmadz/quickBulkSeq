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
#' @param colOrder Character vector, order for columns
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
                         colOrder = NULL,
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

  # For future version: Add params onlySigs = FALSE, onlyNonSigs = FALSE
  # if (onlySigs == TRUE & onlyNonSigs == TRUE) {
  #   stop("Both onlySigs and onlyNonSigs are selected. Please only set one option to TRUE.")
  # }

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
    sigs <- targets
    fontface_vector <- "plain"
  } else {
    sigs <- filtered_res$gene_name[filtered_res$gene_name %in% targets]
  }

  # save for later
  # if (onlySigs) {
  #   targets <- sigs
  # }
  #
  # if (onlyNonSigs) {
  #   targets <- targets[!filtered_res$gene_name %in% targets]
  # }

  numSigs <- nrow(filtered_res)

  selected_genes <- g2name %>%
    filter(gene_name %in% targets) %>%
    arrange(factor(gene_name, levels = targets)) %>%
    distinct(gene_name, .keep_all = TRUE)

  matching_genes <- intersect(selected_genes$gene_id, rownames(assay(vsd)))
  selected_counts <- assay(vsd)[match(matching_genes, rownames(assay(vsd))), ]
  rownames(selected_counts) <- selected_genes$gene_name[match(matching_genes, selected_genes$gene_id)]

  # define z score function
  zscore <- t(scale(t(selected_counts)))

  # determine how significant DEGs are and assign asterisks per level significance
  get_asterisks <- function(gene_name) {
    if (!gene_name %in% filtered_res$gene_name) return("")

    padj <- filtered_res$padj[filtered_res$gene_name == gene_name]

    if (padj < 0.001) return("***")
    if (padj < 0.01) return("**")
    if (padj < 0.05) return("*")
    return("")
  }

  # apply boldface and asterisks to gene labels based on how significant they are
  make_labels <- function(labels, targets) {
    sapply(labels, function(x) {
      asterisks <- get_asterisks(x)
      if (x %in% targets) {
        paste0("bold('", x, asterisks, "')")
      } else {
        paste0("plain('", x, asterisks, "')")
      }
    })
  }

  bold_labels <- make_labels(rownames(zscore), sigs)
  plot_labels <- parse(text = bold_labels)


  if (!is.null(colOrder)) {
    # check if custom column order matches with the actual colnames
    if (!all(colOrder %in% colnames(zscore))) {
      stop(
        paste(
          "Some column names in column_order don't match the data.\nCorrect names are:\n",+paste(colnames(zscore), collapse =
                                                                                                   ", ")
        )
      )
    }
    clustCols <- FALSE  # if custom order, override column clustering
    zscore <- zscore[, colOrder, drop = FALSE]
    annot_info <- annot_info[colOrder, , drop = FALSE]
  }

  if (save) {
    path <- file.path(input$output, input$file_prefix)

    if (!dir.exists(path)) {
      dir.create(path)
    }

    png(
      file.path(path, paste(input$title, "Heatmap.png", sep = " ")),
      width = 2400,
      height = 1800,
      res = 300,
      units = "px"
    )
  }


  final <- pheatmap(
    zscore,
    color = color,
    cluster_rows = clustRows,
    show_rownames = labelRows,
    #labels_row = selected_genes$gene_name,
    labels_row = plot_labels,
    #fontface_row = fontface_vector[match(rownames(zscore)[final$tree_row$order], selected_genes$gene_name)],
    cluster_cols = clustCols,
    border_color = NA,
    annotation_col = annot_info,
    treeheight_row = 0,
    treeheight_col = 0,
    main = paste(input$title, "|", numSigs, "significant genes", sep = " "),  # Title of the heatmap
    silent = TRUE
  )

  print(final)
  if (save) dev.off()

  return (final)
}

