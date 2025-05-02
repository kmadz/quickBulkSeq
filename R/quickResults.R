#' quickResults
#'
#' @param file Either a DESeqDataset object or a path to a DESeqDataset .rds file
#' @param base Character, The type of comparison being made, default "condition"
#' @param num Character, The numerator of the comparison
#' @param denom Character, The denominator of the comparison
#' @param includeNormCounts Logical, whether to include normalized counts in results table, default TRUE
#' @param save Logical, whether to save the file as a .csv to an output path
#' @param title Character, title of the file, defaults to the file name
#' @param output Character, output directory path, default current directory
#' @param filter Logical, whether to filter for significant DEGs
#' @param padj Double, p-adjusted value threshold for filtering if filter = TRUE, default 0.05
#'
#' @importFrom DESeq2 lfcShrink
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr left_join arrange
#' @importFrom magrittr %>%
#' @importFrom readr write_csv
#' @return A results table in dataframe format
#' @export
#'
quickResults <- function(file,
                         base = "condition",
                         title = "",
                         num = "",
                         denom = "",
                         output = "",
                         filter = FALSE,
                         padj = 0.05,
                         save = TRUE,
                         includeNormCounts = TRUE) {
  if (num == "" | denom == "") {
    stop("Numerator and denominator for results not specified")
  }

  input <- ddsHandler(file, output = output, title = title)

  # if (is.character(dds)) {
  #   dds_path <- dds
  #   dds <- readRDS(dds_path)
  #
  #   if (saveToDir) {
  #     output <- dirname(dds_path)
  #     title <- gsub(".*dds_", "", dds_path)
  #     title <- gsub(".rds", "", title)
  #
  #     if (grepl("_", title)) {
  #       title <- gsub("_", " ", title)
  #     }
  #   }
  #
  #   file_prefix <- gsub(pattern = ".rds", replacement = "", file) %>%
  #     gsub(pattern = ".*dds_", replacement = "")
  #
  # } else if (inherits(dds, "DESeqDataSet")) {
  #   file_prefix <- deparse(substitute(file))
  #
  # } else {
  #   stop("file must be either a DESeqDataSet object or a valid file path to an RDS object.")
  # }

  message(paste("Starting results for", input$title, "samples \n", sep = " "))

  # generate results table
  res <- lfcShrink(
    input$dds,
    alpha = 0.05,
    contrast = c(base, num, denom),
    parallel = FALSE,
    type = "ashr"
  )

  # convert results to dataframe and append gene names using g2name
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    left_join(g2name, by = "gene_id") %>%
    arrange(padj)

  if (includeNormCounts) {
    # generate normalized counts dataframe
    norm_df <- as.data.frame(counts(dds, normalized = TRUE)) %>%
      rownames_to_column("gene_id")

    # merge normalized counts to results table
    final_res <- left_join(res_df, norm_df, by = "gene_id")

  } else {
    final_res <- res_df
  }

  # save results as a .csv file
  if (save) {
    path <- file.path(input$output, input$file_prefix)

    if (!dir.exists(path)) {
      dir.create(path)
    }

    write_csv(
      final_res,
      file = file.path(path, paste("RES-", num, "-vs-", denom, ".csv", sep = ""))
    )
  }

  message(paste("Finished results for", input$file_prefix, "samples :)\n", sep = " "))

  if (!filter) return(final_res)

  final_res <- na.omit(final_res)
  final_res <- final_res[final_res$padj < padj, ]

  return (final_res)


}
