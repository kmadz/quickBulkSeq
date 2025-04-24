#' pcaHelper
#'
#' Create a PCA plot that can also adjust for batch effects.
#'
#' @param dds DESeqDataSet object or file path to dds RDS object
#' @param title Character, title of the plot
#' @param label Character, additional text, may be removed
#' @param PC1 Integer, first principal component
#' @param PC2 Integer, second principal component
#' @param intgroup Character or character vector indicating which variables in colData to use for grouping
#' @param batch Logical, apply batch correction using the 'year' column if TRUE
#' @param save Logical, save plot to output directory if TRUE
#' @param output Character, file path to output directory, default is current directory
#' @param saveToDir Logical, extract output directory and title from file path if TRUE
#' @param width Numeric, width of the saved plot in pixels
#' @param height Numeric, height of the saved plot in pixels
#'
#' @importFrom DESeq2 vst plotPCA
#' @importFrom limma removeBatchEffect
#' @importFrom ggplot2 ggsave ggtitle
#' @importFrom ggrepel geom_label_repel
#'
#' @return A PCA ggplot object
#' @export
pcaHelper <- function(dds, title = "", label = "", PC1 = 1, PC2 = 2, intgroup = "condition",
                      blind = FALSE, batch = NULL, save = FALSE, output = ".", saveToDir = FALSE,
                      width = 1600, height = 1200) {

  if (is.character(dds)) {
    dds_path <- dds
    dds <- readRDS(dds_path)
    if (saveToDir) {
      output <- dirname(dds_path)
      title <- gsub(".*dds_", "", dds_path)
      title <- gsub(".rds", "", title)
      if (grepl("_", title)) {
        title <- gsub("_", " ", title)
      }
    }
  } else if (!inherits(dds, "DESeqDataSet")) {
    stop("dds must be either a DESeqDataSet object or a valid file path to an RDS object.")
  }

  vsd <- vst(dds, blind)

  if (!is.null(batch)) {
    assay(vsd) <- removeBatchEffect(assay(vsd), batch = colData(vsd)[[batch]])
  }

  PCA <- plotPCA(vsd, intgroup = intgroup, ntop = 500, pcsToUse = PC1:PC2)
  PCA <- PCA + geom_label_repel(aes(label = name), max.overlaps = 28, force = 3)
  PCA <- PCA + ggtitle(paste(title, " PCA ", PC1, ":", PC2, sep = ""))

  if (save) {
    ggsave(
      filename = file.path(output, paste0(title, label, "PCA.png")),
      plot = PCA,
      width = width,
      height = height,
      dpi = 300,
      units = "px"
    )
  }

  return(PCA)
}
