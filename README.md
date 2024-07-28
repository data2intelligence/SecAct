
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SecAct: Secreted Protein Activity Inference <img src="man/figures/sticker.png" align="right" alt="" width="120" />

<!-- badges: start -->
<!-- badges: end -->

SecAct is an R package desighed for inferring the intercellular activity
of secreted proteins in tumors from gene expression profiles. Users can
input multiple modalities of expression data, including bulk,
single-cell, or spatial transcriptomics. The outputs are 1248 secreted
protein activities for each sample, individual cell, or ST spot,
respectively. SecAct achieves this by leveraging intercellular
signatures trained from 618 genome-wide spatial transcriptomics profiles
across 28 tumor types. If the input are spatial transcriptomics, SecAct
additionally calculates the signaling velocities of secreted proteins
from source to sink spots.

<p align="center">
<img src="man/figures/workflow.png" width="50%"/>
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

## Example

    library(SecAct)

    # prepare expression matrix
    exprPath <- file.path(system.file(package="SecAct"), "extdata/IFNG_GSE100093.diff")
    expr <- read.table(exprPath, sep="\t", check.names=F)

    # run SecAct to infer activity
    res <- SecAct.inference(expr, lambda=10000, nrand=1000)

    # show activity
    head(res$zscore)

## Tutorial

- [Infer the secreted protein activity from gene expression in bulk
  level](https://data2intelligence.github.io/SecAct/articles/Bulk.html)  
- [Infer the secreted protein activity from single-cell RNA-Seq
  data](https://data2intelligence.github.io/SecAct/articles/SingleCell.html)
- [Infer the secreted protein activity from spatial transcriptomics
  data](https://data2intelligence.github.io/SecAct/articles/Spatial.html)

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Kenneth Aldape, Lalage Wakefield,
Peng Jiang. Inferring secreted protein activities at bulk, single-cell,
and spatial levels.
\[<a href="https://github.com/data2intelligence/SecAct" target="_blank">Link</a>\]
