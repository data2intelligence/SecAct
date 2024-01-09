
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SecAct (Secretion Activity) <img src="man/figures/sticker.png" align="right" alt="" width="120" />

<!-- badges: start -->
<!-- badges: end -->

The goal of SecAct is to infer secreted proteins’ activity from RNA-seq
data, including bulk, single-cell, and spatial data.

## Installation

To install `SecAct`, we recommend using `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("data2intelligence/SecAct")
```

Or user can install `SecAct` from the source code. Click
<a href="https://api.github.com/repos/data2intelligence/SecAct/tarball/HEAD" target="_blank">here</a>
to download it.

``` r
# install SecAct in the R environment.
install.packages("Path_to_the_source_code", repos = NULL, type="source")
```

## Example

    library(SecAct)

    Yfile<- file.path(system.file(package = "SecAct"), "extdata/IFNG_GSE100093.diff")
    Y <- read.table(Yfile,sep="\t",check.names=F)

    res <- SecAct.inference(Y, lambda=10000, nrand=1000)

    names(res)
    ## "beta"   "se"     "zscore" "pvalue"

    head(res$beta[,1:2])
    ##           IFNG.15d.AMG811.Lesional.180.mg IFNG.15d.Placebo.Lesional.0.mg
    ## Activin A                   -0.0054571949                   -0.002227756
    ## BDNF                        -0.0165308969                   -0.003590268
    ## BMP2                        -0.0215501262                    0.019876194
    ## BMP4                         0.0039773755                    0.003459698
    ## BMP6                         0.0004974797                    0.025879907
    ## CD40L                       -0.0051185785                    0.001073442

    head(res$zscore[,1:2])
    ##           IFNG.15d.AMG811.Lesional.180.mg IFNG.15d.Placebo.Lesional.0.mg
    ## Activin A                     -1.29207473                     -0.4411385
    ## BDNF                          -3.72453205                     -0.7955953
    ## BMP2                          -4.97028772                      4.5789033
    ## BMP4                           1.04673474                      0.8890111
    ## BMP6                           0.07564311                      5.8074797
    ## CD40L                         -1.25024409                      0.2627249
