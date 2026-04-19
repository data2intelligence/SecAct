
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

To install `SecAct` R package, we recommend using `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("data2intelligence/SecAct")
```

Or user can install `SecAct` R package from the source code. Click
<a href="https://api.github.com/repos/data2intelligence/SecAct/tarball/HEAD" target="_blank">here</a>
to download it.

``` r
# install SecAct dependencies
remotes::install_deps("Path_to_the_source_code", force = TRUE)

# install SecAct R package
install.packages("Path_to_the_source_code", repos = NULL, type="source")
```

The R package has been installed successfully on Operating Systems:

- Red Hat Enterprise Linux 8.10 (Ootpa)
- macOS Sequoia 15.3.1

<div style="border-left: 5px solid #3b82f6; padding: 0px 16px; background-color: #00000; border-radius: 6px;">

##### Alternative installation options:

1.  If you are not familiar with R, you can use the
    <img src="vignettes/img/Python-logo.png" width="2%" style="border:none" />
    Python version of SecAct available
    <a href="https://github.com/data2intelligence/SecActpy" target="_blank">here</a>.

2.  If you prefer not to install SecAct locally, we provide a Docker
    image for both the R and Python versions:

<!-- -->

    # Pull the Docker image.
    docker pull psychemistz/secactpy:with-r

3.  An online server is also available
    <a href="https://appshare.cancer.gov/SecAct/" target="_blank">here</a>.
    Please note that it currently supports only the basic activity
    inference function and is limited to datasets with fewer than 10
    samples.

</div>

## Dependencies

- R version \>= 4.2.0.
- R packages: Matrix, ggplot2, reshape2, patchwork, NMF, akima,
  gganimate, metap, circlize, ComplexHeatmap, ggalluvial, networkD3,
  survival, survminer.
- Optional accelerators (see below):
  [RidgeFast](https://github.com/psychemistz/RidgeFast) (CPU, requires
  GSL) and/or
  [RidgeCuda](https://github.com/psychemistz/RidgeCuda) (GPU, requires
  NVIDIA CUDA Toolkit).
- Optional `rhdf5` for HDF5 streaming of large inputs/outputs in batch
  mode.

## Accelerators and batch mode

The default ridge + permutation kernel is pure R with no compiled
dependencies — SecAct installs and runs anywhere R runs. For large
datasets, two optional accelerator packages provide drop-in speed-ups:

- [RidgeFast](https://github.com/psychemistz/RidgeFast) —
  multi-threaded CPU (GSL + OpenMP)
- [RidgeCuda](https://github.com/psychemistz/RidgeCuda) — NVIDIA GPU
  (CUDA)

With either installed,
`SecAct.activity.inference(..., backend = "auto")` (the default) picks
GPU \> CPU-fast \> pure-R automatically. At `rng_method = "mt19937"`
(default) and `ncores = 1`, output is bit-identical across all three
backends.

For memory-constrained workflows, set `batch_size` to process samples
in column-batches:

``` r
# Memory-bounded inference on 100k samples
res <- SecAct.activity.inference(inputProfile = Y_large,
                                 is.differential = TRUE,
                                 batch_size = 5000)

# Stream results to HDF5 when even the result matrices don't fit in RAM
SecAct.activity.inference(inputProfile = Y_large,
                          is.differential = TRUE,
                          batch_size = 5000,
                          output_h5 = "results.h5",
                          is.group.sig = FALSE)
```

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
activities in intercellular communication. **Nature Methods**, 2026.
\[<a href="https://github.com/data2intelligence/SecAct"
target="_blank">Full Text</a>\]
