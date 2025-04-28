#' quickDESeq
#'
#' Create a DESeqDataSet object using tximport and design terms from a design file.
#'
#' @param design_file CSV file path or data.frame with sample metadata.
#' @param output File path to output directory
#' @param terms Character vector of terms for the model formula.
#' @param mincounts Integer, minimum counts filtering threshold
#' @param filter Integer, minimum count filter threshold.
#' @param save Logical, whether to save the DESeq object as an RDS.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom tximport tximport
#' @importFrom tibble column_to_rownames enframe deframe
#' @importFrom rhdf5 h5closeAll
#' @importFrom DESeq2 DESeqDataSetFromTximport DESeq counts
#' @return A DESeqDataSet object
#' @export
quickDESeq <- function(design_file,
                           pattern = "design_",
                           output  = ".",
                           terms = c("condition"),
                           alignType = "kallisto",
                           save = TRUE,
                           mincounts = 5,
                           filter = 2) {

  # import design file
  sampleData <- tableHandler(design_file, pattern)
  samples <- sampleData$samples
  file_prefix <- sampleData$file_prefix

  # if (csv) {
  #   # gsub searches for the ".csv" file suffix and removes it to extract the file prefix
  #   file_prefix <- gsub(pattern = ".csv", replacement = "", design_file) %>%
  #     gsub(pattern = paste(".*", pattern, sep = ""), replacement = "") #".*design_"
  #   file_prefix
  #   samples <- read.csv(design_file)
  # } else {
  #   # otherwise extract the name from the object presented
  #   name <- deparse(substitute(design_file))
  #   file_prefix <- gsub(pattern = ".*design_", replacement = "", name)
  #   samples <- design_file
  # }

  message(paste("Starting", file_prefix, "samples! \n"))

  if (!"folder" %in% names(samples)) stop("`folder` column not found in design file. Did you add the source folder for your aligned samples?")

  # Add path column for tximport
  samples <- samples %>%
    mutate(path = file.path(folder, file, "abundance.h5"))

  # check if all files specified in path exist
  if (!all(file.exists(samples$path))) {
    missing_files <- samples$path[!file.exists(samples$path)]
    stop(
      "Some files are missing. Make sure the path is correct and the files exist.\n\nFiles missing from path:\n",
      paste(missing_files, collapse = "\n")
    )
  }

  # check if design terms exist
  missing_terms <- setdiff(terms, colnames(samples))
  if (length(missing_terms) > 0) {
    stop(paste("The following design terms are not found in the design file:",
               paste(missing_terms, collapse = ", ")))
  }

  # create tximport object from aligned files
  txi <- tximport(
    files = samples %>% select(sample, path) %>% deframe(),
    type = alignType,
    tx2gene = t2g %>% select(transcript_id, gene_id),
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE
  )

  # rearrange sample data so it's compatible with DESeq
  samples <- samples %>%
    as.data.frame() %>%
    column_to_rownames("sample") %>%
    select(all_of(terms))

  # convert all selected columns to factors
  samples[terms] <- lapply(samples[terms], factor)

  # Build formula
  formula <- as.formula(paste("~", paste(terms, collapse = " + ")))

  # Create DESeqDataSet
  dds <- DESeqDataSetFromTximport(txi, samples, formula)

  # Filtering
  keep <- rowSums(counts(dds) >= mincounts) >= filter
  dds <- dds[keep, ]
  dds <- DESeq(dds)

  # save DESeqDataSet as a .rds file to output directory
  if (save) {
    # check if output directory exists, and if not, create it
    if (!dir.exists(output)) {
      dir.create(output)
    }
    saveRDS(dds, file.path(output, paste0("dds_", file_prefix, ".rds")))
  }

  h5closeAll()
  message("\n done! :)")

  return(dds)
}
