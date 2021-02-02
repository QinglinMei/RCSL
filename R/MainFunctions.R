#'
#' Perform the gene filter
#'
#' @param data the normalized gene expression matrix (each column represents a cell)
#' @param gfRatio the ratio of genes filtering
#'
#' @return gfData gene expression matrix after genes filtering
#' @export
"GenesFilter" <- function(data, gfRatio = 0.025){

  N <- ncol(data)
  G <- nrow(data)
  
  ExpressNumMin <- gfRatio * N
  ExpressNumMax <- (1 - gfRatio) * N
  
  GenesExpress0 <- sapply(1:G, function(n) length(which(data[n,] != 0)))
  gfData <- data[which(GenesExpress0 > ExpressNumMin),,drop = F]
  GenesExpress <- GenesExpress0[which(GenesExpress0 > ExpressNumMin)]
  
  if(max(GenesExpress) == N){
   GenesPre <- gfData[which(GenesExpress > ExpressNumMax),,drop = F]
   GenesVar <- as.matrix(apply(GenesPre, 1, var))
   if(any(GenesVar < 0.9 * mean(GenesVar))){
   gfData <- gfData[-which(GenesVar < 0.9 * mean(GenesVar)),,drop = F]
   }
   
  }
 return(gfData)
}


#'
#' Calculate the initial similarity matrix
#'
#' @param data gene expression matrix after genes filtering (each column represents a cell)
#' @param pcRatio the ratio between the variance of the choosed PCs and the total variance
#' @param largeThre the number of cells to determine whether it is a large datasets
#' @param neiRatio ratio of the number of selected neighbors to the total number of cells
#'
#' @return S initial similarity matrix 
#' @return drData gene expression matrix after PCA processing
#' @export
"SimS" <- function(data, pcRatio = 0.95, gamma = 0.8, largeThre = 2000, NN.method = "LSH", 
                   Dis.method = "Euclidean", LSH.TreeNum = 30, LSH.Dim = 1000, 
                   LSH.Dis = "angular", neiRatio = 0.65){

  res_PCA <- stats::prcomp(t(data), center = TRUE, scale. = TRUE)
  
  if(ncol(data) > 1300){
    EigenValues <- res_PCA$sdev
    V <- sapply(1:length(EigenValues), function(n) sum(EigenValues[1:n])/sum(EigenValues))
    n_dim <- min(which(V >= pcRatio))
    res_pca <- t(res_PCA$x[, 1:n_dim, drop = F])
    if(n_dim>4000){
    res_pca <- t(res_PCA$x[, 1:4000, drop = F])
    }
  }else{
    res_pca <- t(res_PCA$x) 
  }
  
  if(gamma!=1){      
       cat('Calculate the Spearman correlation','\n')
       SpearA <- stats::cor(data, method = "spearman")
       cat('Calculate the Nerighbor Representation','\n')
       NR <- NeigRepresent(res_pca, NN.method = NN.method, Dis.method = Dis.method, 
                        LSH.TreeNum = LSH.TreeNum, LSH.Dim = LSH.Dim, LSH.Dis = LSH.Dis, 
                        neiRatio = neiRatio)
       S <-  gamma * SpearA + (1-gamma) * NR
  } else{
    cat('Calculate the Spearman correlation','\n')
    S <- stats::cor(data, method = "spearman")
  }

  res = list("S" = S, "drData" = res_pca)
  return(res)
}

#'
#' Calculate the neighbor representation of cells 
#' 
#' @import Rcppannoy
#'
#' @param drData gene expression matrix after dimensionality reduced by PCA
#' @param neiRatio ratio of the number of selected neighbors to the total number of cells
#'
#' @return NR the similarity matrix measured by neighbor representation 
#' @export
"NeigRepresent" <- function(drData, NN.method = "KNN", Dis.method = "Euclidean", 
                            LSH.TreeNum = 30, LSH.Dim = 500, LSH.Dis = "angular", 
                            neiRatio = 0.65){

k <- ceiling(ncol(drData) * neiRatio)

if(NN.method == "KNN"){
  if(Dis.method == "Euclidean"){
      cat("Find neighbors by KNN(Euclidean)","\n")
      #distMatrix <- as.matrix(stats::dist(t(drData), diag = T, upper = T))
      distMatrix <- EucDist(drData,drData)
  } 
  if(Dis.method == "Cosine"){
      cat("Find neighbors by KNN(Cosine)","\n")
      Cosine_sim <- t(drData) / sqrt(rowSums(t(drData) * t(drData)))
      Cosine_sim <- Cosine_sim %*% t(Cosine_sim)
      distMatrix <- as.matrix(1 - Cosine_sim)
  }
  neigData <- sapply(1:nrow(distMatrix), function(n){
                     curRow <- sort(distMatrix[n,,drop = T])
                     rownames(as.matrix(which(rank(curRow, ties.method = "first") <= k)))})
}

if(NN.method == "LSH"){
  if(nrow(drData)>LSH.Dim){
    drData <- drData[1:LSH.Dim,,drop = F]
  }
  f = ncol(drData)
  if(LSH.Dis == "angular"){  
      cat("Find neighbors by LSH(Cosine)","\n")
      Tree = methods::new(RcppAnnoy::AnnoyAngular,f)  
  }
  if(LSH.Dis == "hamming"){ 
      cat("Find neighbors by LSH(Hamming)","\n") 
      Tree = methods::new(RcppAnnoy::AnnoyHamming,f)  
  }
  if(LSH.Dis == "Euclidean"){ 
      cat("Find neighbors by LSH(Euclidean)","\n") 
      Tree = methods::new(RcppAnnoy::AnnoyEuclidean,f)  
  }
  for(i in seq(ncol(drData))){
      v = as.vector(drData[,i])
      Tree$addItem(i, v)
  }
  Tree$build(LSH.TreeNum)# 30 trees
  INDEX = mat.or.vec(ncol(drData),k+1)
  for(i in seq(ncol(drData))){
      v = drData[,i]
      INDEX[i,] = Tree$getNNsByVector(v, k+1)
  }
  neigData <- t(INDEX)[-1,]
}

neigData <- neigData[-1,] 
colnames(neigData) <- colnames(drData)

## Solve for reconstruction weights
## Create matrix Z consisting of all neighbors of each cell
wgtsNR <- matrix(0, ncol(drData), ncol(drData))
rownames(wgtsNR) <- colnames(drData)
colnames(wgtsNR) <- colnames(drData)
k <- dim(neigData)[1]

for(i in (1:ncol(neigData))){
  
  # find neighbors of each cell
  neigZ <- drData[, neigData[,i]]
  
  # Calculate the differences between the cell and its neighbors
  diffZ <- neigZ - matrix(drData[,i,drop = F],
                        nrow = nrow(drData),ncol = nrow(neigData),byrow = F)
  localCova <- t(diffZ) %*% diffZ
  
  ## Solve linear system
  if(qr(localCova)$rank < nrow(localCova)){
    wgtsNRi <- pracma::pinv(localCova) %*% matrix(1,nrow=k,ncol=1)
  }else {
    wgtsNRi <- base::solve(localCova, matrix(1,nrow=k,ncol=1))
  }
  
  wgtsNRi <- t(wgtsNRi)
  wgtsNRiRegu <- wgtsNRi / sum(wgtsNRi)
  wgtsNR[i, rownames(localCova)] <- wgtsNRiRegu
 }
 NR <- wgtsNR
 
 return(NR)
}


#'
#' Estimate the optimal number of clusters C for clustering
#' 
#' @param drData gene expression matrix after PCA processing
#' @param S the calculated similarity matrix S from "SimS"
#'
#' @return C the estimated number of clusters
#'
#' @export
"EstClusters" <- function(drData, S){ 

  Diss <- stats::as.dist(1 - S)
  
  if(nrow(S) > 3000){
      
      set.seed(1)
      C <- NbClust::NbClust(t(drData), diss = Diss, distance = NULL,
                            min.nc = 4, max.nc = 12,
                            method = "average", index = "kl")$Best.nc[1]
  }else{ 
      Q0 <- - 10e11
      set.seed(1)
      res <- NbClust::NbClust(t(drData), diss = Diss, distance = NULL,
                              min.nc = 4, max.nc = 12,
                              method = "average", index = "kl")$All.index
      res <- sort(res)
      l <- length(res)
      K_range <- as.numeric(attr(res, "names")[(l - 2):l])
      for(k in K_range){
         resBDSM <- BDSM(S, k)
         a <- resBDSM$y
         one <- matrix(c(1:nrow(S)), byrow = T)
         a <- as.matrix(cbind(one, a))
         b <- matrix(0,nrow = max(a),ncol = max(a))
         for(i in 1:nrow(a)){
             b[a[i,1], a[i,2]] <- 1
         }
         b <- b[, -(which(colSums(b) == 0)),drop = F]
         m <- sum(S) / 2
         c <-  as.matrix(rowSums(S))
         B <- S - (kronecker(matrix(1,1,nrow(c)),c) * kronecker(matrix(1,nrow(c),1),t(c)))/(2 * m)
         Q1 <- 1/(2 * m) * sum(diag(t(b) %*% B %*% b))
         if(Q1 > Q0){
             Q0 <- Q1
             C <- k
         }
       }
	}
  return(C)
}

#' 
#' Solve the problem: min 1/2*x'*L*x-x'*d s.t. x>=0, 1'x=1
#' 
#' @param d 
#' @param l 
#' 
#' @return x 
#'
"EProjSimplexdiag" <- function(d,l){
  lambda <- min(l - d)
  f <- 1
  count <- 1
  
  while(abs(f) > 10^-8){
    v1 <-  1 / l * lambda + d/l
    g <-  sum(1 / l[v1 > 0])
    f <-  sum(v1[v1 > 0]) - 1
    lambda <- lambda - f / g
    if(count > 1000){
    break}
    count <- count + 1
  }
  
  v1 = 1/l*lambda + d/l
  x <- pmax(v1,0)
  
  return(x)
}

#' 
#' Solve the problem: ||A-B||^2 = ||A||^2 + ||B||^2 - 2*A'*B
#' 
#' @param A matrix or vector
#' @param B matrix or vector
#' 
#' @return d matrix or vector
#'
"EucDist" <- function(A, B){
  if (dim(A)[1] == 1){
    A <- rbind(A, matrix(0,1,dim(A)[2]))
    B <- rbind(B, matrix(0,1,dim(B)[2]))
  }
  
  aa <- t(as.matrix(colSums(A * A)))
  bb <- t(as.matrix(colSums(B * B)))
  ab <- t(A) %*% B
  
  d <- kronecker(matrix(1,1,dim(bb)[2]),t(aa)) + kronecker(matrix(1,dim(aa)[2],1), bb) - 2 * ab
  d <-  Re(d)
  d <-  pmax(d,0)
  
  return(d)
}

#' 
#' Calculate the bolock-diagnal matrix B
#' min_{B>=0, B*1=1, F'*F=I}  ||B - A||_1 + r*||B||^2 + 2*lambda*trace(F'*L*F)
#' 
#' @import igraph
#'
#' @param S the calculated initial similarity matrix S
#' @param C the estimated number of clusters C
#' 
#' @return B block-diagonal matrix
#' @return y clustering results
#' @export
"BDSM" <- function(S, C){
  NITER <- 30
  zr  <-  10e-11
  lambda  <-  0.1
  r  <-  0

  S  <-  S - diag(diag(S))
  num <- nrow(S)
  A10 <-  (S + t(S)) / 2
  L0  <-  diag(colSums(A10)) - A10
  L0 <- pmax(L0, t(L0))
  tmp0 <- eigen(L0) 
  w0 <- tmp0$values
  V0 <- tmp0$vectors
  wsort0 <- sort(w0, index.return = TRUE)
  di0 <- wsort0$x
  idx0 <- wsort0$ix
  idx1 <- idx0[1 : num]
  F  <-  V0[, idx1][, 1:C]
 
  if(sum(w0[1 : (C + 1)]) < zr){
     stop('The original graph has more than %d connected component', C)
   }
  
  if(sum(w0[1 : C]) < zr){
     g0 <- igraph::graph_from_adjacency_matrix(A10)
     y <- as.matrix(components(g0, mode = "strong")$membership, drop = F)
     B <-  S
     return
   }
    
  u = vector('list', num)  
  
  for(j in 1:num){
     s0 <-  S[j,,drop = F]
     idxa0 <- 1:num
     u[[j]] <-  matrix(1, 1, length(idxa0))
  }
 
  for(iter in 1:NITER){
      dist0 <- EucDist(t(F), t(F))
      B <-  matrix(0, num, num)
      for(i in 1:num){
         idxa0 <-  1:num
         ai <-  S[i,,drop = F][idxa0]
         di <-  dist0[i,idxa0,drop = F]
         ad <- u[[i]]*ai-lambda*di/2
         si <- EProjSimplexdiag(ad, u[[i]] + r * rep(1,length(idxa0)))
         u[[i]] <-  1 / (2 * sqrt((si - ai)^2 + .Machine$double.eps))
         B[i,idxa0] <-  si 
        }
      A <-  B
      A <-  (A + t(A)) / 2
      L <-  diag(colSums(A)) - A
      F_old <-  F
      L <-  pmax(L, t(L))
      w <-  eigen(L)$values
      V <-  eigen(L)$vectors
      wsort <-  sort(w,index.return = TRUE)
      idx2 <- wsort$ix[1:C]
      F <-  V[, idx2]
      ev <-  w[wsort$ix]
      fn1 <- sum(ev[1:C])
      fn2 <- sum(ev[1:C + 1])
      if(fn1 > zr){ lambda <- 2 * lambda}
        else if(fn2 < zr) {
                   lambda <- lambda / 2 
                   F <- F_old
        }
       else {break}
   }
   cat("======== Calculate maximal strongly connected components ========","\n")
   g <- igraph::graph_from_adjacency_matrix(B,mode = "undirected",weighted = T,diag = FALSE)
   G <- igraph::components(g, mode = "strong"); 
   y <- as.matrix(G$membership,drop = F)
   clusternum <- G$no
   if(clusternum != C){
     sprintf('Can not find the correct cluster number: %d', C)
    }
   res <- list("B" = B, "y" = y)
   return(res)
}


   
   