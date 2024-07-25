
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SecAct: Secreted Protein Activity Inference <img src="man/figures/sticker.png" align="right" alt="" width="120" />

<!-- badges: start -->
<!-- badges: end -->

SecAct is an R package desighed for inferring the intercellular activity
of secreted proteins in tumors from gene expression profiles. Users can
input multiple modalities of expression data, including bulk,
single-cell, or spatial transcriptomics. The outputs are 1248 secreted
protein activities for each sample, individual cell, or ST spot,
respectively. If input data are spatial transcriptomics, SecAct
additionally calculates the secreted signaling velocity at each spatial
spot. The velocity direction starts from the source cell producing a
secreted protein and moves to sink cells receiving the secreted protein
signal. The velocity magnitude represents the product between the
secreted protein-coding gene expression at source cells and signaling
activities at sink cells.

<p align="center">
<img src="man/figures/workflow.png" width="52%"/>
</p>

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

## Dependencies

- R version \>= 4.2.0.
- R packages: Matrix, ggplot2, patchwork.
- C Library: GSL.

## Example 1 (Input: expression matrix)

    library(SecAct)

    Yfile <- file.path(system.file(package = "SecAct"), "extdata/IFNG_GSE100093.diff")
    Y <- read.table(Yfile,sep="\t",check.names=F)

    res <- SecAct.inference(Y, lambda=10000, nrand=1000)

    head(res$zscore[,1:2])
    ##           IFNG.15d.AMG811.Lesional.180.mg IFNG.15d.Placebo.Lesional.0.mg
    ## Activin A                      -1.2489955                     -0.5632170
    ## BDNF                           -3.7953628                     -0.7830596
    ## BMP2                           -4.9343815                      4.4209498
    ## BMP4                            0.9777525                      0.8402344
    ## BMP6                            0.1167711                      5.7487538
    ## CD40L                          -1.2868588                      0.3233309

## Tutorial

- [Bulk](https://data2intelligence.github.io/SecAct/articles/Bulk.html)  
- [Single-Cell](https://data2intelligence.github.io/SecAct/articles/SingleCell.html)
- [Spatial](https://data2intelligence.github.io/SecAct/articles/Spatial.html)
