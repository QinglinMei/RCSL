# ## Construct simulated datasets 

# library(splatter)
# library(scater)

# ## ============= 300 cells 1 =============
# NCells <- 300
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 10000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.25,0.25,0.3), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# data <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 300 cells 2 =============
# NCells <- 300
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 15000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.25,0.2,0.2,0.15), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# data <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 500 cells 1 =============
# NCells <- 500
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 10000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.2,0.3,0.3), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 500 cells 2 =============
# NCells <- 500
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 15000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.3,0.2,0.2,0.1), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 1000 cells 1 =============
# NCells <- 1000
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 10000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.2,0.3,0.2,0.1), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 1000 cells 2 =============
# NCells <- 1000
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 15000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.15,0.1,0.25,0.2,0.1,0.2), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 2000 cells 1 =============
# NCells <- 2000
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 10000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.2,0.3,0.2,0.1), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 2000 cells 2 =============
# NCells <- 2000
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 15000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.1,0.2,0.2,0.1,0.2), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 3000 cells 1 =============
# NCells <- 3000
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 10000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.2,0.2,0.15,0.2,0.1,0.15), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

# ## ============= 3000 cells 2 =============
# NCells <- 3000
# DataName <- paste0('SimulateData','_',NCells)
# set.seed(1)
# sce <- mockSCE(ncells = NCells, ngenes = 15000)
# params <- splatEstimate(sce)
# sim.groups <- splatSimulate(params, group.prob=c(0.15,0.1,0.25,0.15,0.1,0.15,0.1), method = "groups",
#                             verbose = FALSE)
# sim.groups <- logNormCounts(sim.groups)
# curData <- logcounts(sim.groups)
# label <-  sim.groups$Group

