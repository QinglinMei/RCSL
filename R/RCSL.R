#'
#' Perform the RCSL program
#'
#' @param datat normalizaed gene expression matrix(each column represents a cell)
#' @param GF should I need the gene filter step?
#' @param gfRatio the ratio of the gene filter
#' @param pcRatio the ratio between the variance of the choosed PCs and the total variance
#' @param largeThre the number of cells to determine whether it is a large datasets
#' @param neiRatio ratio of the number of selected neighbors to the total number of cells in NR
#'
#' @return gfData gene expression matrix after genes filtering
#' @return B block-diagonal matrix
#' @return C estimated number of clusters
#' @return y clustering results
#'

RCSL <- function(data, GF = TRUE, gfRatio = 0.025, pcRatio = 0.95, largeThre = 2000,
                 neiRatio = 0.65){

# Gene filter
if(GF == TRUE){
gfData <- GenesFilter(data, gfRatio)
} else{
gfData <- data[rowSums(data)!=0,,drop=F]
}

# Calculate the similarity matrix S
resSimS <- SimS(gfData, pcRatio, largeThre)

# Estimate the number of clusters C
Estimated_C <- EstClusters(resSimS$drData,resSimS$S)

# Calculate the block diagonal matrix B
res_BDSM <- BDSM(resSimS$S, Estimated_C)

B <- res_BDSM$B
y <- res_BDSM$y

res <- list("gfData" = gfData, "drData" = resSimS$drData,
            "S" = resSimS$S, "B" = B, "C" = C, "y" = y)
return(res)
}





