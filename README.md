# Rank Constrained Similarity Learning (RCSL)
RCSL is an R toolkit for single-cell clustering and trajectory analysis using single-cell RNA-seq data.

## Installation
This package can be insatlled through devtools in R:
```{r}
$ R
> library("devtools")
> devtools::install_github("QinglinMei/RCSL")
```
Now RCSL can be loaded in R:
```{r}
> library(RCSL)
```
## Input

The input to RCSL is a normalized data matrix with columns being cells and rows being genes in log(CPM+1), log(RPKM+1), log(TPM+1) or log(FPKM+1) format. 

## Usage

We provide an example script to run RCSL in *demo_RCSL.R*. 

There are also nine functions available in RCSL serving as sub-commands.

Subcommand | Description
-----------|----------
`GenesFilter` | Perform genes filtering.
`SimS` | Calculate the initial similarity matrix S.
`NeigRepresent` | Calculate the neighbor representation of cells.
`EstClusters` | Estimate the optimal number of clusters C.
`BDSM` | Learn the block-diognal matrix B.
`PlotMST` | Construct MST based on clustering results from RCSL.
`PlotPseudoTime` | Infer the pseudo-temporal ordering of cells.
`getLineage` | Infer the lineage based on the clustering results and the starting cell.
`PlotTrajectory` | Plot the developmental trajectory based on the clustering results and the starting cell.

## Required annotations for RCSL

1) The full RCSL workflow requires extra tools: our tool requires three R packages to run, namely the *igraph* package (see https://igraph.org/) to find the stronggest connected components, the *ggplot2* package (see https://cran.r-project.org/web/packages/ggplot2/index.html) to plot the developmental trajectory and MST, the *SingleCellExperiment* package (see https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) to read our *SingleCellExperiment* object.
2) Five example datasets provided in the directory *Data* were downloaded from https://hemberg-lab.github.io/scRNA.seq.datasets/. These datasets are quality controlled by scater criterion and stored in SingleCellExperiment Bioconductor S4 classes.


## DEBUG

Please feel free to contact us if you have problems running our tool at meiql@mail.sdu.edu.cn.




