#' DESeqHelperTxi
#'
#' Create a DESeqDataSet object using tximport and design terms from a design file.
#'
#' @param design_file CSV file path or data.frame with sample metadata.
#' @param output File path to output directory
#' @param terms Character vector of terms for the model formula.
#' @param csv Logical, TRUE if design_file is a CSV.
#' @param mincounts Integer, minimum counts filtering threshold
#' @param filter Integer, minimum count filter threshold.
#' @param save Logical, whether to save the DESeq object as an RDS.
#'
#' @return A DESeqDataSet object
#' @export
DESeqHelperTxi <- function(design_file, output  = ".", terms = c("condition"), csv = TRUE, save = TRUE, mincounts = 5, filter = 2) {

  # import design file
  if (csv) {
    # gsub searches for the ".csv" file suffix and removes it to extract the file prefix
    file_prefix <- gsub(pattern = ".csv", replacement = "", design_file) %>%
      gsub(pattern = ".*design_", replacement = "")
    samples <- read.csv(file.path(main, design_file))
  } else {
    # otherwise extract the name from the object presented
    name <- deparse(substitute(design_file))
    file_prefix <- gsub(pattern = ".*design", replacement = "", name)
    samples <- design_file
  }

  message(paste("Starting", file_prefix, "samples! \n"))

  # Add path column for tximport
  samples <- samples %>%
    mutate(path = file.path(folder, file, "abundance.h5"))

  # Validate selected design terms
  missing_terms <- setdiff(design_terms, colnames(samples))
  if (length(missing_terms) > 0) {
    stop(paste("The following design terms are not found in the design file:",
               paste(missing_terms, collapse = ", ")))
  }

  # Prepare tximport
  txi <- tximport(
    files = samples %>% select(sample, path) %>% deframe(),
    type = "kallisto",
    tx2gene = t2g %>% select(transcript_id, gene_id),
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE
  )

  # Subset sample data for DESeq
  samples <- samples %>%
    as.data.frame() %>%
    column_to_rownames("sample") %>%
    select(all_of(design_terms))

  # Convert all selected columns to factors
  samples[design_terms] <- lapply(samples[design_terms], factor)

  # Build formula
  formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))

  # Create DESeqDataSet
  dds <- DESeqDataSetFromTximport(txi, samples, formula)

  # Filtering
  keep <- rowSums(counts(dds) >= mincounts) >= filter
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  # Save object
  if (!dir.exists(output)) {
    dir.create(output)
  }

  if (save) {
    saveRDS(dds, file.path(output, paste0("dds_", file_prefix, ".rds")))
  }

  h5closeAll()
  message("done! :) \n")

  return(dds)
}
