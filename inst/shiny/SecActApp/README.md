# SecAct Unified Visualization App

Interactive Shiny application for exploring secreted protein activity across
bulk, single-cell, and spatial transcriptomics data.

## Quick Start

```r
# Install SecAct (if not already installed)
devtools::install_github("data2intelligence/SecAct")

# Launch the app
SecAct::runSecActApp()
```

## Tabs

- **Spatial** — Visualize SecAct activity on spatial coordinates (Visium, CosMx, Xenium, etc.)
- Bulk, Single Cell, and Inference tabs coming in future releases.

## Requirements

- R >= 3.5.0
- SpaCET package (for spatial visualization)
- shiny, shinyWidgets, bslib, DT (installed automatically via Suggests)
- Optional: SeuratObject, SeuratDisk (for .h5seurat file support)
