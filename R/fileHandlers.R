#' ddsHandler
#'
#' @param dds
#'
#' @return a list specific to the type of object
#' @export
#'
ddsHandler <- function(file, output = ".") {
  if (is.null(file) | is.na(file)) {
    stop("Please enter a filepath to a DESeqDataSet .rds file OR a DESeqDataSet object, input was either null or NA")
  }


  if (is.character(dds)) {
    dds <- readRDS(file)

    if (output == ".") {
      output <- dirname(dds_path)
    }

    title <- gsub(".*dds_", "", dds_path)
    title <- gsub(".rds", "", title)

    file_prefix <- title

    if (grepl("_", title)) {
      title <- gsub("_", " ", title)
    }

    # return a list with information exclusive to this instance of the function
    final <- list(
      dds = dds,
      output = output,
      title = title,
      file_prefix = file_prefix
    )

    return(final)

  } else if (inherits(dds, "DESeqDataSet")) {
    dds <- file

    # since dds is a variable name, use deparse(substitute(x)) to extract the name
    title <- gsub(".*dds_", "", deparse(substitute(file)))

    file_prefix <- title

    if (grepl("_", title)) {
      title <- gsub("_", " ", title)
    }

    final <- list(
      dds = dds,
      output = output,
      title = title,
      file_prefix = file_prefix
    )

    return(final)
  } else {
    stop("file must be either a DESeqDataSet object or a valid file path to an RDS object.")
  }
}

#' tableHandler
#'
#' @param design_file a design file of either .csv or dataframe type
#'
#' @return a list with a samples table and a file prefix
#' @export
#'
tableHandler <- function(design_file, pattern = "design_") {
  if (is.character(design_file)) {

    file_prefix <- gsub(pattern = ".csv", replacement = "", design_file) %>%
      gsub(pattern = paste(".*", pattern, sep = ""),
           replacement = "") #".*design_"

    samples <- read.csv(design_file)

    final <- list(samples = samples, file_prefix = file_prefix)
    return(final)

  } else {

    # otherwise extract the name from the object presented
    name <- deparse(substitute(design_file))
    file_prefix <- gsub(pattern = paste(".*", pattern, sep = ""),
                        replacement = "",
                        name)

    samples <- design_file

    final <- list(samples = samples, file_prefix = file_prefix)

    return (final)
  }
}
