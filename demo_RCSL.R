
setwd("./RCSL/Data")

library(SingleCellExperiment)
library(ggplot2)


## Load dataset 
# For example, the input data "Goolam.rds"

origData <- readRDS("./Data/Goolam.rds")
data <- logcounts(origData)
label <- origData$cell_type1
DataName <- "Goolam"

# In this example, the input data is could be also "Goolam.txt"

## Run RCSL
# RCSL in Detail
# The main `RCSL` method is a wrapper that calls other `RCSL` functions in the following order:
# GenesFilter(); SimS(); EstClusters(); BDSM()

library(RCSL)

res_RCSL <- RCSL(data, GF = TRUE, gfRatio = 0.025, pcRatio = 0.95, largeThre = 2000)



## Trajectory Analysis to time-series datasets
# TrajectoryAnalysis in Detail
# The main `TrajectoryAnalysis` method also calls other functions in the following order:
# PlotPseudoTime(); getLineage(); PlotMST(); PlotTrajectory()
# The MST plot can be displayed without the starting cell. 
# The pseudo-temporal ordering and the developmental trajectories need the starting cell as input.

res_TrajecAnalysis <- TrajectoryAnalysis(res_RCSL$gfData, res_RCSL$drData, res_RCSL$S,
                                         clustRes = res_RCSL$y, TrueLabel = label, startPoint = 1,
                                         dataName = DataName)
										 



# Display the constructed MST 

res_TrajecAnalysis$MSTPlot


# Display the plot of the pseudo-temporal ordering 

res_TrajecAnalysis$PseudoTimePlot


# Display the plot of the inferred developmental trajectory

res_TrajecAnalysis$TrajectoryPlot


# The clustering results can be evaluted by ARI

ARI_RCSL <- igraph::compare(res_RCSL$y, label, method = "adjusted.rand")
