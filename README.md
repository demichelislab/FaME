# FaME

FAME is a tool that can test the aberrations of the samples in one or more
datasets in order to find signals of mutual exclusivity or co-occurrence
between all the pairs of genes. FaME can test many different combinations of
aberrations very quickly by using several techniques described in the methods
section of the paper.

## Running FaME

> FaME requires large amount of computational resources if run
> on large datasets.

In order to run FaME just copy the genomic data that can be downloaded from
[here](https://github.com/demichelislab/SPICE-pipeline#pre-computed-data)
to the data/genomic/ folder and run the code. The script requires the R
*librarian* package to be installed. The *librarian* package can be installed
with the following command:

```r
install.packages('remotes')
remotes::install_github('https://github.com/DesiQuintans/librarian')
```

In order to run the analysis just run the script `fame_analysis.R` with the
following command:

```bash
Rscript code/fame_analysis.R
```

When launched the script will automatically install the packages that are
needed and then the analysis will start.

## OpenBLAS

FaME leverage R matrix multiplication. The code to compute this operation is
provided by the BLAS library. Even though the default matrix multiplication
works out of the box in order to speed up the computation is better to install
the OpenBLAS library. Such library provide a much faster version of the
operation. If you have an ubuntu system you can follow
[this guide](https://www.r-bloggers.com/2013/07/for-faster-r-use-openblas-instead-better-than-atlas-trivial-to-switch-to-on-ubuntu/)
to install OpenBLAS or to check that the library is installed.

