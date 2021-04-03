#'
#' Plot the visualization of constructed Minimum Spanning Tree based on the clustering results of RCSL
#'
#' @param drData preprocessed gene expression data 
#‘               (each column represent a cell)
#' @param clustRes the clustering results identified by RCSL
#' @param TrueLabel the real cell types to color the dots in plot
#' @param dataName the name of the data that will be showed in the plot
#' @param fontSize the font size of the plot
#' @param VisualMethod the method for 2D visualization including UMAP,t-SNE and PCA
#'
#' @return MSTPlot ggplot object of the visualization of constructed MST
#'
#' @export
#'
#' @examples
#' gfData <- GenesFilter(yan)
#' TrueLabel <- ann$cell_type1
#' res_SimS <- SimS(gfData)
#' C <- EstClusters(res_SimS$drData,res_SimS$S)
#' res_BDSM <- BDSM(res_SimS$S,C)
#' PlotMST(res_SimS$drData,res_BDSM$y)
#'
"PlotMST" <- function(drData, clustRes, TrueLabel, dataName = "",
                    fontSize = 12, VisualMethod = "umap"){

  centers_num <- max(clustRes)
  centers <- matrix(nrow=nrow(drData),ncol = centers_num)

  for(i in seq_len(centers_num)){
      centers[,i] <- rowMeans(drData[,which(clustRes == i),drop = FALSE])
  }

  centersSim <- (stats::cor(centers, method = "kendall") + 1) / 2
  g <- igraph::graph_from_adjacency_matrix(centersSim,mode = "undirected",weighted = TRUE,diag = FALSE)
  g_mst <- igraph::mst(g)

  if(VisualMethod == "umap"){
     reducedData <- data.frame(umap::umap(t(drData),n_components = 2)$layout)
     x.lab = "umap1"
     y.lab = "umap2"
  }
  if(VisualMethod == "tsne"){
     reducedData <- data.frame(Rtsne::Rtsne(t(drData),perplexity = 10)$Y)
     x.lab = "tsne1"
     y.lab = "tsne2"
  }
  if(VisualMethod == "pca"){
     res_PCA <- stats::prcomp(t(drData),center = TRUE, scale. = TRUE)$x
     reducedData <- data.frame(res_PCA[,seq_len(2)])
     colnames(reducedData) <- c("X1","X2")
     x.lab = "pca1"
     y.lab = "pca2"
  }
  reducedData.centers <- matrix(nrow = max(clustRes),ncol = 2)
  for(i in seq_len(max(clustRes))){
     reducedData.centers[i,] = colMeans(reducedData[which(clustRes == i),])
  }
  gg <- igraph::get.data.frame(g_mst)
  df_centers <- as.data.frame(reducedData.centers)
  df_centers$vertices <- c(seq_len(nrow(reducedData.centers)))
  gg$from.x <- df_centers$V1[match(gg$from, df_centers$vertices)]
  gg$from.y <- df_centers$V2[match(gg$from, df_centers$vertices)]
  gg$to.x <- df_centers$V1[match(gg$to, df_centers$vertices)]
  gg$to.y <- df_centers$V2[match(gg$to, df_centers$vertices)]
  reducedData$label <- TrueLabel
  colorFeature <- c("#f54291","#1089ff","#f3c623","#7e0cf5","#fb832d","#32ff6a",
                    "#B07AA1","#B65412","#FF6666","#7c3c21","#ff427f","#56C188")
  MSTPlot <- ggplot2::ggplot() + geom_point(data = reducedData, aes(x = X1, y = X2, color = label),size = 1.3,alpha = 0.75) +
    scale_colour_manual(values = colorFeature) +
    geom_segment(data = gg,aes(x = from.x, xend = to.x, y = from.y, yend = to.y),
                 size = 1.2,colour = "black") +
    geom_point(data = df_centers,aes(x = V1, y = V2),size = 2.3,alpha = 0.8,colour = "black") +
    ggtitle(dataName) + labs(x = x.lab, y = y.lab) +
    theme(panel.background = element_blank(),
          legend.text = element_text(size = (fontSize - 2)),
          axis.title = element_text(size = fontSize),
          axis.text = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="right",legend.title = element_blank(),
          plot.title = element_text(size = fontSize, face = "bold", hjust = 0.5),
          legend.key = element_blank(), legend.background = element_blank(),
          legend.box.margin=margin(-15,0,0,-9),
          legend.justification = "top",
          legend.key.width = unit(0, "lines"),
          legend.key.size = unit(1,"lines"),
          plot.margin = unit(c(0.5, 1.5, 1.5, 0.5), "lines"))+
          guides(colour = guide_legend(direction = "vertical"))

  return(MSTPlot)
}

#'
#' Infer the pseudo-temporal ordering between the cell types using the
#' distance from a cell type to the predefined starting cell type.
#'
#' @param S the similarity matrix calculated by SimS() function
#' @param TrueLabel the real cell types used to indicate the vertical axis
#' @param startPoint the posiition of the starting cell in the matrix
#' @param fontSize the font size of the plot
#' @param dataName the name of the data that will be showed in the plot
#' @param sim indicate the input data is simialrity matrix or not
#'
#' @return PstudoTime
#' @return PseudoTimePlot  ggplot object of the pseudo-temporal ordering of cells
#'
#' @export
#'
#' @examples
#' gfData <- GenesFilter(yan)
#' TrueLabel <- ann$cell_type1
#' res_SimS <- SimS(gfData)
#' C <- EstClusters(res_SimS$drData,res_SimS$S)
#' res_BDSM <- BDSM(res_SimS$S,C)
#' PlotPseudoTime(res_SimS$S,res_BDSM$C,staetpoint=1)
#'
"PlotPseudoTime" <- function(S, TrueLabel, startPoint, fontSize = 12,
                             dataName = "", sim = TRUE){

  if(sim == TRUE){
     NormA <- S / sum(S)
     Dist <- (1 - NormA)
  } else{
     Dist = as.matrix(stats::dist(S))
  }

  startPoint <- startPoint+1

  pseudoTime <- Dist[startPoint,]

  gInput <- data.frame( x = pseudoTime[-1],
                        y = TrueLabel[-1],
                        Timepoint = TrueLabel[-1])

  colorFeature <- c("#f54291","#1089ff","#f3c623","#7e0cf5","#fb832d","#32ff6a",
                    "#B07AA1","#B65412","#FF6666","#7c3c21","#ff427f","#56C188")

  PseudoTimePlot <- ggplot2::ggplot(gInput, aes(x = x, y = y, colour = Timepoint)) +
       geom_point(size = 1.1,alpha = 0.78) + labs(x = "Pseudo Time") +
       scale_colour_manual(values = colorFeature) +
       ggtitle(dataName)+
       theme( panel.background = element_blank(), panel.border = element_blank(),
              panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = fontSize,colour = "black"),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size = fontSize),
              plot.title = element_text(size = fontSize, face = "bold", hjust = 0.5),
              legend.text = element_text(size = (fontSize)),
              legend.position = "none",legend.title = element_blank(),legend.box = "vertical",
              legend.key = element_blank(),
              legend.spacing.x = unit(0, 'cm'))

  res <- list("PstudoTime" = pseudoTime, "PseudoTimePlot" = PseudoTimePlot )

  return(res)
}

#'
#' Infer the development lineage based on the clustering results from RCSL and the pseudotime
#'
#' @param drData preprocessed gene expression data (each column represent a cell)
#' @param clustRes the clustering results identified by RCSL
#' @param pseudoTime inferred by PlotPseudoTime() using the similarity matrix S and starting cell
#' @param simMeasure the calculation method of measuring the cluster centers' similarity
#'
#' @return lineage the cell lineages connected all the cluster centers based on the clustering results from RCSL
#'
#' @export
#'
#' @examples
#' gfData <- GenesFilter(yan)
#' TrueLabel <- ann$cell_type1
#' res_SimS <- SimS(gfData)
#' C <- EstClusters(res_SimS$drData,res_SimS$S)
#' res_BDSM <- BDSM(res_SimS$S,C)
#' Pseudo <- PlotPseudoTime(res_SimS$S,res_BDSM$C,staetpoint=1)
#' getLineage(res_SimS$drData,res_BDSM$y,Pseudo$pseudoTime)
#'
"getLineage" <- function(drData, clustRes, pseudoTime, simMeasure = "kendall"){
  centers_num <- max(clustRes)
  centers <- matrix(nrow = nrow(drData),ncol = centers_num)
  centersPseudo <- numeric(length(centers_num))
  for(i in seq_len(centers_num)){
     centers[,i] <- rowMeans(drData[,which(clustRes == i),drop = FALSE])
     centersPseudo[i] <- mean(pseudoTime[which(clustRes == i)])
  }

  if(simMeasure == "spearman"){
     centersSim <- (stats::cor(centers, method = "spearman") + 1) / 2
  }

  if(simMeasure == "pearson"){
     centersSim <- (stats::cor(centers, method = "pearson") + 1) / 2
  }

  if(simMeasure == "kendall"){
     centersSim <- (stats::cor(centers, method = "kendall") + 1) / 2
  }

  for(j in seq_len(max(clustRes))){
     centersPseudo[j] <- mean(pseudoTime[which(clustRes == j),drop = FALSE])
  }
  g <- igraph::graph_from_adjacency_matrix(centersSim,mode = "undirected",weighted = TRUE,diag = FALSE)
  g_mst <- igraph::mst(g)
  lineage <- numeric()
  Edges <- igraph::get.edgelist(g_mst)
  lineagePre <- table(Edges)
  LineageNum <- max(lineagePre) - 1

  if(LineageNum == 1){
     lineage[c(1,centers_num)] <- which(lineagePre == 1)
     AdjacentNodes = as.list(igraph::adjacent_vertices(g_mst,v = c(seq_len(centers_num))))
     for(i in seq_len(centers_num - 2)){
         AdjNodes = as.numeric(AdjacentNodes[[lineage[i]]])
         lineage[i + 1] = AdjNodes[which(!(AdjNodes %in% lineage))]
     }
  }else{
     vertices <- c(1:centers_num)
     longestPath <- as.numeric(get_diameter(g_mst))
     vers <- vertices[which(!vertices %in% longestPath)]
     EdgesRest <- matrix(nrow = length(vers), ncol = 2)
     for(x in seq_len(length(vers))){
         EdgesRest[x,] <- Edges[which(Edges == vers[x], arr.ind = TRUE)[1],]
      }
     nei <- EdgesRest[EdgesRest != vers]
     for(i in seq_len(length(nei))){
         neis <- longestPath[c(which(longestPath == nei[i]) - 1,which(longestPath == nei[i]) + 1)]
         nearest <- neis[which.max(centersSim[neis, vers[i]])]
         if(which(longestPath == nei[i]) < which(longestPath == nearest)){
             lineage <- append(longestPath, vers, after = which(longestPath == nei[i]))
         }else{
             lineage <- append(longestPath, vers, after = (which(longestPath == nei[i]) - 1))
         }
     }
    lineage = rank(centersPseudo[lineage])
  }

  return(lineage)
}

#'
#' Infer the developmental trajectories based on the clustering results from RCSL
#'
#' @param gfData preprocessed gene expression data (each column represent a cell)
#' @param clustRes the clustering results identified by RCSL
#' @param TrueLabel the real cell types
#' @param lineage  the lineage obtained by getLineage()
#' @param fontSize the size of font in the plot
#' @param dataName the name of the data that will be showed in the plot
#' @param VisualMethod the display method of 2-D visualization
#'
#' @return TrajectoryPlot ggplot object of the inferred developmental trajectories
#'
#' @export
#'
#' @examples
#' gfData <- GenesFilter(yan)
#' TrueLabel <- ann$cell_type1
#' res_SimS <- SimS(gfData)
#' C <- EstClusters(res_SimS$drData,res_SimS$S)
#' res_BDSM <- BDSM(res_SimS$S,C)
#' Linea <- getLineage(res_SimS$drData,res_BDSM$y,Pseudo$pseudoTime)
#' PlotTrajectory(gfData,res_BDSM$y,TrueLabel,lineage=Linea)
#'
"PlotTrajectory" <- function(gfData, clustRes, TrueLabel, lineage,
                             fontSize = 12, dataName = "", VisualMethod = "umap"){
  if(VisualMethod == "umap"){
     reducedData <- data.frame(umap::umap(t(gfData),n_components = 2)$layout)
     x.lab = "Umap1"
     y.lab = "Umap2"
  }

  if(VisualMethod == "tsne"){
     reducedData <- data.frame(Rtsne::Rtsne(t(gfData),perplexity = 10)$Y)
     x.lab = "tsne1"
     y.lab = "tsne2"
  }

  if(VisualMethod == "pca"){
     res_PCA <- stats::prcomp(t(gfData),center = TRUE, scale. = TRUE)$x
     reducedData <- data.frame(res_PCA[,seq_len(2)])
     colnames(reducedData) <- c("X1","X2")
     x.lab = "Pca1"
     y.lab = "Pca2"
  }

  reducedDataCenters <- matrix(nrow = max(clustRes), ncol = 2)
  for(i in seq_len(max(clustRes))){
     reducedDataCenters[i,] = colMeans(reducedData[which(clustRes == i),])
  }

  reducedData$label <- TrueLabel
  reducedDataCenters <- as.data.frame(reducedDataCenters[lineage,,drop = FALSE])
  plot(0)
  ps <- data.frame(xspline(reducedDataCenters, shape = -0.2, lwd = 2, draw = FALSE))
  dev.off()

  colorFeature <- c("#f54291","#1089ff","#f3c623","#7e0cf5","#fb832d","#32ff6a",
                    "#B07AA1","#B65412","#FF6666","#7c3c21","#ff427f","#56C188")

  TrajectoryPlot <- ggplot2::ggplot() + geom_point(data = reducedData, aes(x = X1, y = X2, color = label),size = 1.1,alpha = 0.75) +
       scale_colour_manual(values = colorFeature) +
       geom_path(data = ps, aes(x, y),  colour = "black",size = 1.3,
                 arrow = arrow(length = unit(0.3, "cm"), angle = 30, type = "closed"),
                 lineend = "round",linejoin = "round") +
       geom_smooth()+
       ggtitle(dataName)+
       guides(color = guide_legend(nrow = 1)) +
       labs(x = x.lab, y = y.lab) +
       theme(axis.text = element_blank(),
             panel.background = element_blank(), panel.border = element_blank(),
             panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
             axis.line = element_line(colour = "black"),
             axis.title = element_text(size = fontSize),
             plot.title = element_text(size = fontSize, face = "bold", hjust = 0.5),
             legend.position = "right",legend.title = element_blank(),
             legend.text = element_text(size = (fontSize - 2)),
             legend.key = element_blank(), legend.background = element_blank(),
             legend.box.margin=margin(-15,0,0,-9),
             legend.justification = "top",
             legend.key.width = unit(0, "lines"),
             legend.key.size = unit(1,"lines"))+
			 guides(colour = guide_legend(direction = "vertical"))
  return(TrajectoryPlot)
}

#'
#' Trajectory analysis 
#'
#' @param gfData preprocessed gene expression data (each column represent a cell)
#' @param drData preprocessed gene expression data (each column represent a cell)
#' @param S the similarity matrix calculated by SimS() function
#' @param clustRes the clustering results identified by RCSL
#' @param fontSize the size of font in the plot
#' @param TrueLabel the real cell types used to indicate the vertical axis
#' @param startPoint the posiition of the starting cell in the matrix
#' @param dataName the name of the data that will be showed in the plot
#' @param simMeasure the calculation method of measuring the cluster centers' similarity
#' @param VisualMethod the display method of 2-D visualization
#'
#' @return PseudoTimePlot, MSTPlot, TrajectoryPlot
#'
#' @export
#'
#' @examples
#' gfData <- GenesFilter(yan)
#' TrueLabel <- ann$cell_type1
#' res_SimS <- SimS(gfData)
#' C <- EstClusters(res_SimS$drData,res_SimS$S)
#' res_BDSM <- BDSM(res_SimS$S,C)
#' TrajectoryAnalysis(gfData,res_SimS$drData,res_SimS$S,res_BDSM$y,
#'                    TrueLabel=TrueLabel,startPoint=1)
#'
"TrajectoryAnalysis" <- function(gfData, drData, S, clustRes, fontSize = 12, TrueLabel,
                                 startPoint, dataName = "", sim = TRUE,
                                 simMeasure = "kendall", VisualMethod = "umap"){

  PseudoRes <- PlotPseudoTime(S, TrueLabel, startPoint, fontSize, dataName)

  lineage <- getLineage(drData, clustRes, PseudoRes$PstudoTime, simMeasure = "kendall")

  PseudoTimePlot <- PseudoRes$PseudoTimePlot

  MSTPlot <- PlotMST(drData, clustRes, TrueLabel, dataName, fontSize, VisualMethod)

  TrajectoryPlot <- PlotTrajectory(gfData, clustRes, TrueLabel, lineage, fontSize,
                                    dataName, VisualMethod)
  res <- list("PseudoTimePlot" = PseudoTimePlot,
              "MSTPlot" = MSTPlot,
              "TrajectoryPlot" = TrajectoryPlot)

  return(res)
}
