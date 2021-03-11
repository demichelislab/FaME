# FaME

FAME is a tool that can test the aberrations of the samples in one or more
datasets in order to find signals of mutual exclusivity or co-occurrence
between all the pairs of genes. FaME can test many different combinations of
aberrations very quickly by using several techniques described in the methods
section of the paper.

## Running FaME

> :warning: FaME requires large amount of computational resources if run
> on large datasets.

In order to run FaME just copy the genomic data that can be downloaded from the
links reported above to the data/genomic/ folder and run the code. The script
requires the R *librarian* package. The *librarian* package can be installed with
the following command:

```r
install.packages('librarian')
```

In order to run the analysis just run the script `fame_analysis.R`. When
launched the script will automatically install the packages that are needed and
then the analysis will start.
