#' quickVolcano
#'
#' @param results a table of results of either a tibble format or csv
#' @param title Character, the name of the plot to be shown
#' @param save Logical, whether to save the plot
#' @param FCcutoff Double, cutoff for log2FoldChange threshold, default 1
#' @param pCutoff Double, cutoff for p-adjusted value for significance, default 0.05
#' @param output Character, path to save image to, default current directory
#' @param targets Character vector, a list of target genes to display instead of the top n genes, default NULL
#' @param xlim Integer vector in form of c(lowerbound, upperbound), defines x limit for volcano plot
#' @param ylim Integer vector in form of c(lowerbound, upperbound), defines x limit for volcano plot
#' @param numTopGenes Number of top genes to display
#' @param labelSize Double, size of labels
#' @param displayNumSigs Logical, whether to show how many significant genes, default TRUE
#'
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom ggplot2 ggsave
#'
#' @return A volcano plot
#' @export
#'
quickVolcano <- function(results,
                         title = "",
                         output = "",
                         targets = "",
                         xlim = c(min(results[["log2FoldChange"]], na.rm = TRUE) - 1.5, max(results[["log2FoldChange"]], na.rm = TRUE) +
                                    1.5),
                         ylim = c(0, max(-log10(results[["padj"]]), na.rm = TRUE) + 5),
                         numTopGenes = 10,
                         labelSize = 4,
                         save = FALSE,
                         displayNumSigs = TRUE,
                         FCcutoff = 1,
                         pCutoff = 0.05) {
  input <- resultsHandler(results, title = title, output = output)

  res <- input$results

  filtered_res <- na.omit(res)
  filtered_res <- filtered_res[filtered_res$padj < pCutoff, ]

  if (displayNumSigs) {
    numUp <- nrow(filtered_res[filtered_res$log2FoldChange > 0, ])
    numDown <- nrow(filtered_res[filtered_res$log2FoldChange < 0, ])
  }

  if (targets == "") {
    top_genes <- filtered_res[order(filtered_res$padj), ][1:numTopGenes, ]
    top_genes <- top_genes$gene_id
  } else {
    top_genes <- targets
  }

  final <- EnhancedVolcano(
    res,
    lab = ifelse(res$gene_id %in% top_genes, res$gene_name, ""),
    x = 'log2FoldChange',
    y = 'padj',
    xlim = xlim,
    ylim = ylim,
    ylab = 'FDR',
    pCutoff = pCutoff,
    FCcutoff = FCcutoff,
    pointSize = 1.5,
    labSize = labelSize,
    title = title,
    subtitle = ifelse(
      displayNumSigs == TRUE,
      paste(
        numUp,
        'upregulated DEGs,',
        numDown,
        'downregulated DEGs',
        sep = " "
      ),
      ""
    ),
    caption = paste('FC cutoff = ', FCcutoff, '; FDR cutoff = ', pCutoff),
    legendPosition = "right",
    legendLabSize = 14,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.9,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    max.overlaps = Inf,
    boxedLabels = TRUE
  )

  plot(final)

  if (save) {
    path <- file.path(input$output, input$file_prefix)

    if (!dir.exists(path)) {
      dir.create(path)
    }

    ggsave(
      file.path(path, paste(input$title, "Volcano.png", sep = " ")),
      width = 2400,
      height = 1800,
      dpi = 300,
      units = "px"
    )
  }
  return(final)
}
