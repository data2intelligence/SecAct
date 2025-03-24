
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SecAct: Secreted Protein Activity Inference <img src="man/figures/sticker.png" align="right" alt="" width="120" />

<!-- badges: start -->
<!-- badges: end -->

SecAct is an R package designed for inferring the intercellular
signaling activity of secreted proteins from gene expression profiles.
Users can input multiple modalities of expression data, including
spatial, single-cell, or bulk transcriptomics data. The outputs are the
inferred signaling activities of 1248 secreted proteins for each spatial
spot, individual cell, or sample, depending on the input data type.
Based on the inferred activities, SecAct provides multiple downstream
application modules. For spatial data, SecAct can infer the signaling
pattern and signaling velocity for secreted proteins. For single-cell
data, SecAct can infer the intercellular communication network and
signaling flow from source cells to receiver cells. For bulk data,
SecAct can infer secreted protein risk scores for a large cohort linked
to clinical data, and can infer secreted protein activities that are
differentially regulated between two phenotypes. These functionalities
and terms are explained more formally in the following tutorials.

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
# install SecAct in the R environment.
install.packages("Path_to_the_source_code", repos = NULL, type="source")
```

The package has been installed successfully on Operating Systems:

- Red Hat Enterprise Linux 8.10 (Ootpa)
- macOS Sequoia 15.3.1

<details>
<summary>
How to install GNU Scientific Library (GSL)?
</summary>

If your operating system does not already have GSL installed, please
follow one of the installation methods below, depending on your
operating system.

    # Linux (Ubuntu/Debian):
    sudo apt-get install libgsl-dev

    # Linux (Fedora/RHEL-based):
    sudo dnf install gsl-devel

    # macOS (Homebrew):
    brew install gsl

</details>

## Dependencies

- C Library: GNU Scientific Library (GSL).
- R version \>= 4.2.0.
- R packages: Matrix, ggplot2, reshape2, patchwork, NMF, gganimate,
  metap, circlize, ComplexHeatmap, ggalluvial, networkD3, survival,
  survminer.

## Example

``` r
library(SecAct)

dataPath <- file.path(system.file(package = "SecAct"), "extdata/")
expr.diff <- read.table(paste0(dataPath, "Ly86-Fc_vs_Vehicle_logFC.txt"))
res <- SecAct.activity.inference(inputProfile=expr.diff, is.differential=TRUE)

head(res$zscore)

##            Change
## A1BG     7.588670
## A2M     15.974063
## A2ML1    3.303167
## AADACL2 -1.789684
## ABHD15  -5.353521
## ABI3BP  18.785484
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

- [Cell-cell communication mediated by secreted
  proteins](https://data2intelligence.github.io/SecAct/articles/scCCC.html)  
- [Secreted protein signaling activity for distinct cell
  states](https://data2intelligence.github.io/SecAct/articles/scState.html)

#### Bulk RNA sequencing data

- [Clinical relevance of secreted proteins in a large patient
  cohort](https://data2intelligence.github.io/SecAct/articles/bulkCohort.html)
- [Secreted protein signaling activity change between two
  phenotypes](https://data2intelligence.github.io/SecAct/articles/bulkChange.html)

## Citation

Beibei Ru, Lanqi Gong, Emily Yang, Kenneth Aldape, Lalage Wakefield,
Peng Jiang. Inference of secreted protein activities in intercellular
communication.
\[<a href="https://github.com/data2intelligence/SecAct" target="_blank">Link</a>\]
