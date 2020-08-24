# Rank Constrained Similarity Learning (RCSL)
RCSL is an R toolkit for single cell clustering and trajectory analysis using RNA-seq data.

## Installation
This package can be insatlled through devtools in R:
```{r}
$ R
> library("devtools")
> devtools::install_github("QinglinMei/RCSL")
```
Now we can load RCSL:
```{r}
> library(RCSL)
```
## Input

The input of RCSL is expected that the preprocessing steps like sequencing quality control (QC) or normalisation are performed by a user in advance. Then, the input is assumed to be an log(CPM+1) (or log(RPKM+1), log(TPM+1), log(FPKM+1)) matrix with columns correspond to cells and genes correspond to genes. 

## Usage

We provide an *R* demo code to run *RCSL* in the script *demo_RCSL.R*. 

There are also nine functions available in *RCSL* serving as sub-commands.

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

### Required annotations for RCSL

1) The full RCSL workflow requires extra tools: our tool requires three R packages to run, namely the *igraph* package (see https://igraph.org/) to find the stronggest connected components, the *ggplot2* package (see https://cran.r-project.org/web/packages/ggplot2/index.html) to plot the developmental trajectory and MST, the *SingleCellExperiment* package (see https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) to read our *SingleCellExperiment* object.
2) Five example datasets provided in the directory *Data* were downloaded from https://hemberg-lab.github.io/scRNA.seq.datasets/. These datasets are quality controlled by scater criterion and stored in SingleCellExperiment Bioconductor S4 classes.


### DEBUG

Please feel free to contact us if you have problems running our tool at meiql@mail.sdu.edu.cn.




