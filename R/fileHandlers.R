#' tableHandler
#'
#' @param design_file a design file of either .csv or dataframe type
#'
#' @return a list with a samples table and a file prefix
#' @export
#'
designHandler <- function(design_file, pattern = "design_") {
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

#' ddsHandler
#'
#' @param file Filepath to a DESeqDataSet .rds file OR a DESeqDataSet object
#' @param output Path to output directory, default current directory
#'
#' @return a list specific to the type of object
#' @export
#'
ddsHandler <- function(file, output = ".", title = "") {
  # if (is.null(file) | is.na(file)) {
  #   stop("Please enter a filepath to a DESeqDataSet .rds file OR a DESeqDataSet object, input was either null or NA")
  # }

  if (is.character(file)) {
    dds <- readRDS(file)

    if (output == "") {
      output <- dirname(file)
    }

    if (title == "") {
      title <- gsub(".*dds_", "", file)
      title <- gsub(".rds", "", title)
    }

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

  } else if (inherits(file, "DESeqDataSet")) {
    dds <- file

    if (title == "") {
      # since dds is a variable name, use deparse(substitute(x)) to extract the name
      title <- gsub(".*dds_", "", deparse(substitute(file)))
    }

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

#' resultsHandler
#'
#' @param results a path to a results table in the form of "RES-name.csv" or a results table dataframe
#' @param title title of the results table, defaults to deparsed name
#' @param output path to save what this table is used for, defaults to current directory
#'
#' @return A list containing a results table, title, and output
#' @export
#'
resultsHandler <- function(results,
                           title = "",
                           output = ".") {
  if (is.data.frame(results)) {
    if (title == "") {
      title <- deparse(substitute(results))
      title <- gsub("_", " ", title)
    }
  } else if (is.character(results)) {
    if (output == "") {
      output <- dirname(results)
    }
    results <- as.data.frame(read.csv(results))
    if (title == "") {
      title <- gsub(".*RES-", "", results) %>%
        gsub(pattern = ".csv", replace = "") %>%
        gsub(pattern = "-", replace = " ") %>%
        gsub(pattern = "_", replace = " ")
    }
  }

  final <- list(results = results,
                title = title,
                output = output)

}

