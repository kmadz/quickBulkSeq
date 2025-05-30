% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quickResults.R
\name{quickResults}
\alias{quickResults}
\title{quickResults}
\usage{
quickResults(
  file,
  base = "condition",
  title = "",
  num = "",
  denom = "",
  output = "",
  filter = FALSE,
  padj = 0.05,
  save = TRUE,
  includeNormCounts = TRUE
)
}
\arguments{
\item{file}{Either a DESeqDataset object or a path to a DESeqDataset .rds file}

\item{base}{Character, The type of comparison being made, default "condition"}

\item{title}{Character, title of the file, defaults to the file name}

\item{num}{Character, The numerator of the comparison}

\item{denom}{Character, The denominator of the comparison}

\item{output}{Character, output directory path, default current directory}

\item{filter}{Logical, whether to filter for significant DEGs}

\item{padj}{Double, p-adjusted value threshold for filtering if filter = TRUE, default 0.05}

\item{save}{Logical, whether to save the file as a .csv to an output path}

\item{includeNormCounts}{Logical, whether to include normalized counts in results table, default TRUE}
}
\value{
A results table in dataframe format
}
\description{
quickResults
}
