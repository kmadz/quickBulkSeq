pcaHelper <- function(dds, title = "",  label = "", PC1 = 1, PC2 = 2, save = TRUE, output = ".", saveToDir = FALSE) {
  if (saveToDir) {
    if (class(dds) == "character") {
      dir <- dirname(dds)
      title <- gsub(".*dds_", "", dds)
      title <- gsub(".rds", "", title)

      if (grepl("_", title)) {
        title <- gsub("_", " ", title)
      }
    } else {
      stop("Error in saveToDir: dds is not a character and cannot have its directory extracted")
    }
  } else {
    dir <- output
  }

  dds <- readRDS(dds)
  vsd <- vst(dds, blind = FALSE)

  assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch = vsd$year)

  z <- plotPCA(vsd,intgroup=c("condition"),ntop=500, pcsToUse = PC1:PC2)
  z <- z + geom_label_repel(aes(label = name), max.overlaps = 28, force = 3)
  z <- z + ggtitle(label = paste(title, " PCA ", PC1, ":", PC2, sep = ""))
  # plot(z)
  #
  # ggsave(
  #   paste(dir, "/", title, label, " PCA.png", sep = ""),
  #   width = width,
  #   height = height,
  #   dpi = 300,
  #   units = "px"
  # )
  return(z)
}
