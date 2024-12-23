
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SecAct: Secreted Protein Activity Inference <img src="man/figures/sticker.png" align="right" alt="" width="120" />

<!-- badges: start -->
<!-- badges: end -->

SecAct is an R package desighed for inferring the intercellular activity
of secreted proteins in tumors from gene expression profiles. Users can
input multiple modalities of expression data, including bulk,
single-cell, or spatial transcriptomics (ST). The outputs are 1248
secreted protein activities for each sample, individual cell, or ST
spot, respectively. SecAct achieves this by leveraging intercellular
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

The R package has been installed successfully on Operating systems:

- macOS Sonoma 14.5
- Rocky Linux 8.7 (Green Obsidian)

## Example

    library(SecAct)

    # prepare expression matrix
    exprPath <- file.path(system.file(package="SecAct"), "extdata/IFNG_GSE100093.diff.gz")
    expr <- read.csv(exprPath, row.names=1, check.names=F)

    # run SecAct to infer activity; ~3 mins
    res <- SecAct.inference(expr)

    # show activity
    head(res$zscore)

    ##         Anti-IFNG.15day Anti-IFNG.57day
    ## A1BG        27.41700606       15.560484
    ## A2M         -1.00132544        3.703776
    ## A2ML1       -5.70139021       -5.585635
    ## AADACL2     -0.03412573        8.556714
    ## ABHD15       7.91357240        3.360326
    ## ABI3BP      -1.57170759       -4.292832

    # show IFNG activity
    res$zscore["IFNG",]

    ## Anti-IFNG.15day Anti-IFNG.57day 
    ##       -10.76148       -30.50151 

## Tutorial

#### Spatial transcriptomcis (ST) data

- [Signaling patterns and velocities for multi-cellular ST
  data](https://data2intelligence.github.io/SecAct/articles/stPattern.html)  
- [Intercellular communication for single-cell resolution ST
  data](https://data2intelligence.github.io/SecAct/articles/stCCC.html)

#### Single-cell RNA sequencing data

- [Cell-cell communication mediated by secreted
  proteins](https://data2intelligence.github.io/SecAct/articles/scCCC.html)  
- [Secreted protein signaling activity for various cell
  states](https://data2intelligence.github.io/SecAct/articles/scState.html)

#### Bulk RNA sequencing data

- [Secreted protein activity difference between two
  phenotypes](https://data2intelligence.github.io/SecAct/articles/bulkDiff.html)
- [Clinical relavance of secreted proteins in a large patient
  cohort](https://data2intelligence.github.io/SecAct/articles/bulkCohort.html)

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Kenneth Aldape, Lalage Wakefield,
Peng Jiang. Inference of secreted protein activities in intercellular
communication.
\[<a href="https://github.com/data2intelligence/SecAct" target="_blank">Link</a>\]
