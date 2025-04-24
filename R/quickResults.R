#' quickResults
#'
#' @param file Either a DESeqDataset object or a path to a DESeqDataset .rds file
#' @param base Character, The type of comparison being made, default "condition"
#' @param num Character, The numerator of the comparison
#' @param denom Character, The denominator of the comparison
#' @param includeNormCounts Logical, whether to include normalized counts in results table, default TRUE
#' @param save Logical, whether to save the file as a .csv to an output path
#'
#' @importFrom DESeq2 lfcShrink
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr left_join arrange
#' @importFrom magrittr %>%
#' @return A results table in dataframe format
#' @export
#'
quickResults <- function(file,
                         base = "condition",
                         num = "",
                         denom = "",
                         output = ".",
                         save = TRUE,
                         includeNormCounts = TRUE) {
  if (num == "" | denom == "") {
    stop("Numerator and denominator for results not specified")
  }

  file_prefix <- gsub(pattern = ".rds", replacement = "", file) %>%
    gsub(pattern = ".*dds_", replacement = "")

  path <- file.path(output, file_prefix)

  if (!dir.exists(path)) {
    dir.create(path)
  }

  message(paste("Starting results for", file_prefix, "samples \n", sep = " "))

  dds <- readRDS(file)

  # generate results table
  res <- lfcShrink(
    dds,
    alpha = 0.05,
    contrast = c(base, num, denom),
    parallel = FALSE,
    type = "ashr"
  )

  # convert results to dataframe and append gene names using t2g
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    left_join(g2name, by = "gene_id") %>%
    arrange(pvalue)

  if (includeNormCounts) {
    # generate normalized counts dataframe
    norm_df <- as.data.frame(counts(dds, normalized = TRUE)) %>%
      rownames_to_column("gene_id")

    # merge normalized counts to results table
    final_res <- left_join(res_df, norm_df, by = "gene_id")

  } else {
    final_res <- res_df
  }

  if (!dir.exists(path)) {
    dir.create(file.path(path))
  }

  # save results as a .csv file
  if (save) {
    write_csv(
      final_res,
      file = paste(
        path,
        "/",
        file_prefix,
        ".",
        num,
        ".vs.",
        denom,
        ".RESULTS.csv",
        sep = ""
      )
    )
  }

  message(paste("Finished results for", file_prefix, "samples :)\n", sep = " "))

  return(final_res)
}
