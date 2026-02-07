
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SecAct: Secreted Protein Activity Inference <img src="man/figures/sticker.png" align="right" alt="" width="120" />

<!-- badges: start -->

<!-- badges: end -->

SecAct is an R package designed for inferring the intercellular
signaling activity of secreted proteins from gene expression profiles.
Users can input multiple modalities of expression data, including
spatial, single-cell, or bulk transcriptomics data. The outputs are the
inferred signaling activities of \>1,000 secreted proteins for each
spatial spot, individual cell, or sample, depending on the input data
type. Based on the inferred activities, SecAct provides multiple
downstream application modules. For spatial data, SecAct can infer the
signaling pattern and signaling velocity for secreted proteins. For
single-cell data, SecAct can infer the intercellular communication
network and signaling flow from source cells to receiver cells. For bulk
data, SecAct can infer secreted protein risk scores for a large cohort
linked to clinical data, and can infer secreted protein activities that
are differentially regulated between two phenotypes. These
functionalities and terms are explained more formally in the following
tutorials.

<p align="center">

<img src="man/figures/workflow.png" width="100%"/>
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
# install.packages("remotes")
remotes::install_deps("Path_to_the_source_code", force = TRUE)

# install SecAct in the R environment.
install.packages("Path_to_the_source_code", repos = NULL, type="source")
```

The package has been installed successfully on Operating Systems:

- Red Hat Enterprise Linux 8.10 (Ootpa)
- macOS Sequoia 15.3.1

## Dependencies

- C Library: GNU Scientific Library (GSL).
- R version \>= 4.2.0.
- R packages: Matrix, ggplot2, reshape2, patchwork, NMF, akima,
  gganimate, metap, circlize, ComplexHeatmap, ggalluvial, networkD3,
  survival, survminer.

## Example

``` r
library(SecAct)

dataPath <- file.path(system.file(package = "SecAct"), "extdata/")
expr.diff <- read.table(paste0(dataPath, "Ly86-Fc_vs_Vehicle_logFC.txt"))

# infer activity; ~2 mins
res <- SecAct.activity.inference(inputProfile=expr.diff, is.differential=TRUE) 

head(res$zscore)
```

## Tutorial

SecAct is applicable to multiple modalities of gene expression profiles,
including spatial, single-cell, and bulk transcriptomics data. The
following tutorials demonstrate its applications across each data type.

#### Spatial transcriptomcis (ST) data

- [Signaling patterns and velocities for multi-cellular ST
  data](https://data2intelligence.github.io/SecAct/articles/stPattern.html)  
- [Intercellular communication for single-cell resolution ST
  data](https://data2intelligence.github.io/SecAct/articles/stCCC.html)

#### Single-cell RNA sequencing data

- [Secreted protein signaling activity for distinct cell
  states](https://data2intelligence.github.io/SecAct/articles/scState.html)
- [Cell-cell communication mediated by secreted
  proteins](https://data2intelligence.github.io/SecAct/articles/scCCC.html)

#### Bulk RNA sequencing data

- [Secreted protein signaling activity change between two
  phenotypes](https://data2intelligence.github.io/SecAct/articles/bulkChange.html)
- [Clinical relevance of secreted proteins in a large patient
  cohort](https://data2intelligence.github.io/SecAct/articles/bulkCohort.html)

## Contact

For questions, bug reports, or feature requests, please submit an
[issue](https://github.com/data2intelligence/SecAct/issues). To keep the
issue tracker focused and constructive, advertising or promotional
content is not permitted.

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Seongyong Park, George Zaki, Kenneth
Aldape, Lalage Wakefield, Peng Jiang. Inference of secreted protein
activities in intercellular communication.
\[<a href="https://github.com/data2intelligence/SecAct" target="_blank">Link</a>\]
