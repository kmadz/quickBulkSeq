# quickBulkSeq

quickBulkSeq is an R package designed to simplify the process of analyzing bulk RNA-sequencing data. The package comes with the ability to generate differential expression tables, volcano plots, and heatmaps. Currently, quickBulkSeq uses DESeq and tximport to import and analyze transcript-level abundances of RNA-seq reads, but more features such as importing gene-level abundances using featureCounts will come in the future.

## Installation

Run the following code to install quickBulkSeq.

```R
library(devtools)
install_github("kmadz/quickBulkSeq")
```

## Usage

```R
library(qBulkSeq)

# returns a DESeqDataSet
quickDESeq(design_table_or_file, alignType)

# returns a DESeq results table
quickResults(dds_object_or_file, numerator, denominator)

# returns a volcano plot, can also highlight specific genes if desired
quickVolcano(results_table_or_file, title, targets)

# return a heatmap, can also select only specified genes
quickHeatmap(dds_object_or_file, numerator, denominator, title, targets)
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
