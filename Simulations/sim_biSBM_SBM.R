# Script to compare biSBM to SBM

rm(list=ls())

library(ggplot2)
library(cowplot)
library(pheatmap)
library(dplyr)
library(fossil)

source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")


getFDR <- function(inferredGraph,binaryTruth){
  truthVec <- as.vector(binaryTruth)
  inferredGraphVec <- as.vector(inferredGraph)
  FDR <- sum(inferredGraphVec[truthVec==0])/ max(sum(inferredGraphVec),1)
  TDR <- sum(inferredGraphVec[truthVec==1])/ sum(truthVec==1)
  res <- list(FDR=FDR, TDR=TDR)
  return(res)
}

OutFolder <- "../Output/biSBM/"

mu <- as.numeric(commandArgs(trailingOnly=T)[1])
jid <- as.numeric(commandArgs(trailingOnly=T)[2])

cat('mean parameter is: ', mu, '. \n')
today <- "20250611"
n1 <- 40
n2 <- 60
Q1 <- 2
Q2 <- 3
network <- 'bisbm'
model <- 'Gauss01'

ICL <- 0
exclusive_Rows_Cols <- F #ifelse(Q1==Q2, T, F)

# generate modular bipartite network ----
nu <- matrix(mu,Q1,Q2)  #%>% `diag<-`(., 1)
w <- matrix(0.1,Q1,Q2) %>% `diag<-`(., 0.8)

if (network == 'sbm'){
  n <- n1+n2
  Q <- Q1+Q2
  theta <- list(pi= rep(1/Q,Q), nu0=c(0,1))
  nu <- rbind(cbind(matrix(mu, Q1, Q1), nu), 
              cbind(t(nu), matrix(mu, Q2, Q2)))
  theta$nu <- matrix(c(nu[lower.tri(nu,diag = T)],
                       rep(0.5,Q*(Q+1)/2)),
                     nrow=Q*(Q+1)/2,ncol=2)
  w <- rbind(cbind(matrix(0, Q1, Q1), w), 
             cbind(t(w), matrix(0, Q2, Q2)))
  theta$w <- w[lower.tri(w,diag = TRUE)]
  
  obs <- noisySBM::rnsbm(n, theta, modelFamily='Gauss')

  A.full <- obs$latentAdj[order(obs$latentZ), order(obs$latentZ)]
  row_clusters <- obs$latentZ[order(obs$latentZ)]
  n1 <- sum(row_clusters %in% 1:Q1)
  n2 <- sum(row_clusters %in% (Q1+1):(Q1+Q2))
  # pheatmap(A.full, cluster_rows = F, cluster_cols = F)
  A <- A.full[(row_clusters %in% 1:Q1), (row_clusters %in% (Q1+1):(Q1+Q2))]
  pheatmap(A, cluster_rows = F, cluster_cols = F)
  R <- obs$dataMatrix[order(obs$latentZ), order(obs$latentZ)]
  # pheatmap(R, cluster_rows = F, cluster_cols = F)
  dataMatrix <- R[(row_clusters %in% 1:Q1), (row_clusters %in% (Q1+1):(Q1+Q2))]
  # pheatmap(dataMatrix, cluster_rows = F, cluster_cols = F)
}

if (network == 'bisbm'){
  theta <- list(alpha1 = rep(1/Q1,Q1),
                alpha2 = rep(1/Q2,Q2),
                nu0 = c(0,1))
  theta$nu <- list(mean = nu,
                   sd = matrix(0.5,Q1,Q2))
  theta$pi <- w

  # theta$pi
  obs <- rnbisbm(n1,n2,theta = theta)
  A <- obs$latentAdj
  A.full <- rbind(cbind(matrix(0,n1,n1),A),
                  cbind(t(A), matrix(0,n2,n2)))
  
  dataMatrix <- obs$dataMatrix
  R <- rbind(cbind(matrix(rnorm(n1^2),n1,n1),obs$dataMatrix),
             cbind(t(obs$dataMatrix), matrix(rnorm(n2^2),n2,n2)))
  
  # pheatmap(A[order(obs$latentZ[[1]]), order(obs$latentZ[[2]])], cluster_rows = F, cluster_cols = F)
}

# gene

# generate preferential attachment network ----
if (network == 'pa'){
  # A <- readRDS(paste0(OutFolder,'network_',paste0(today,"_n1_",n1,"_n2_",n2,"_Q1_",Q1,"_Q2_",Q2, '_mu3'), '.rds'))
  A <- readRDS(paste0(OutFolder,'network_pa_n1_26n2_75.rds'))
  n1 <- nrow(A)
  
  dataMatrix <- matrix(rnorm(n1*n2), n1, n2)
  dataMatrix[A==1] <- rnorm(sum(A==1), mean = mu)
  A.full <- rbind(cbind(matrix(0,n1,n1),A),
                  cbind(t(A), matrix(0,n2,n2)))
  
  R <- rbind(cbind(matrix(rnorm(n1^2),n1,n1),dataMatrix),
             cbind(t(dataMatrix), matrix(rnorm(n2^2),n2,n2)))
}
file.end <- paste0(today,"_n_",n1+n2,"_Q1_",Q1,"_Q2_",Q2, '_mu', mu, '_ICL', ICL)

file.end <- paste0(file.end, '_', network)

dataVector <- as.vector(dataMatrix)
pvals <- 2*pnorm(abs(dataVector),lower.tail = FALSE)

# run noisy biSBM ----
start.time = Sys.time()
res.bisbm <- fitNobiSBM(dataMatrix, model=model,
                        exclusive_Rows_Cols = F, 
                        sbmSize = if(ICL) list(Q1=1:5, Q2=1:5, explor=1.5) else list(Q1=Q1, Q2=Q2, explor=1.5),
                        initParam = list(nbOfbeta=1, nbOfPointsPerbeta=NULL,
                                         maxNbOfPasses=2, minNbOfPasses=1),
                        nbCores = 1)
end.time = Sys.time()
bisbm.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

# plotICL(res.bisbm, exclusive_Rows_Cols)
best.id <- which(sapply(res.bisbm,function(a) a$sbmParam$ICL)==getBestQ(res.bisbm, exclusive_Rows_Cols)$ICL)
currentSolution.bisbm <- res.bisbm[[best.id]]

resGraph.bisbm <- graphInferenceNobiSBM(dataMatrix,
                                        currentSolution.bisbm$clustering_row,
                                        currentSolution.bisbm$clustering_col,
                                        currentSolution.bisbm$theta, alpha=0.1)

cat('noisybiSBM', unlist(getFDR(resGraph.bisbm$A, A)), '...\n')


alpha.range <- c(0.005,0.025,0.05,0.1,0.15,0.25)
FDRTDR <- list()

FDRTDR[['randIndexbiSBM']] <- list('row'=fossil::adj.rand.index(currentSolution.bisbm$clustering_row, obs$latentZ[[1]]),
                                   'col'=fossil::adj.rand.index(currentSolution.bisbm$clustering_col, obs$latentZ[[2]]))

FDRTDR[['biSBM']] <- sapply(alpha.range, function(a) {
  Ahat <- matrix(1*(resGraph.bisbm$qvalues<a), n1, n2)
  FDR <- sum(Ahat[A==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A==1])/ sum(A==1)
  data.frame(FDR=FDR,TDR=TDR)
})
colnames(FDRTDR[['biSBM']]) <- alpha.range

# run noisySBM ----
library(noisySBM)
start.time = Sys.time()
res.sbm <- fitNSBM(R, model = model,
                   sbmSize = if(ICL) list(Qmin = 1, Qmax = Q1+Q2+1, explor = 1.5) else  list(Qmin = Q1+Q2, Qmax = Q1+Q2, explor = 1.5),
                   initParam = list(nbOfTau = 1, nbOfPointsPerTau = NULL,
                                    maxNbOfPasses = 2,
                                    minNbOfPasses = 1))
end.time = Sys.time()
sbm.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

# noisySBM::plotICL(res.sbm)
best.id <- which(sapply(res.sbm,function(a) a$sbmParam$ICL)==noisySBM::getBestQ(res.sbm)$ICL)
currentSolution.sbm <- res.sbm[[best.id]]

resGraph.sbm <- noisySBM::graphInference(R,
                                         currentSolution.sbm$clustering,
                                         currentSolution.sbm$theta, alpha=0.1)

cat('noisySBM', unlist(getFDR(resGraph.sbm$A, A.full)), '...\n')

# pheatmap(resGraph.sbm$A, cluster_rows = F, cluster_cols = F)

FDRTDR[['SBM']] <- sapply(alpha.range, function(a) {
  Ahat <- matrix(0, n1+n2, n1+n2)
  Ahat[lower.tri(Ahat)] <- (resGraph.sbm$qvalues < a)
  Ahat <- Ahat + t(Ahat)
  FDR <- sum(Ahat[A.full==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A.full==1])/ sum(A.full==1)
  data.frame(FDR=FDR,TDR=TDR)
})
colnames(FDRTDR[['SBM']]) <- alpha.range

FDRTDR[['randIndexSBM']] <- list('row'=fossil::adj.rand.index(currentSolution.sbm$clustering[1:n1], obs$latentZ[[1]]),
                                   'col'=fossil::adj.rand.index(currentSolution.sbm$clustering[-(1:n1)], obs$latentZ[[2]]))


FDRTDR[['SBM.2']] <- sapply(alpha.range, function(a) {
  Ahat <- matrix(0, n1+n2, n1+n2)
  Ahat[lower.tri(Ahat)] <- (resGraph.sbm$qvalues < a)
  Ahat <- Ahat + t(Ahat)
  Ahat <- Ahat[1:n1, -(1:n1)]
  FDR <- sum(Ahat[A==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A==1])/ sum(A==1)
  data.frame(FDR=FDR,TDR=TDR)
})
colnames(FDRTDR[['SBM.2']]) <- alpha.range

FDRTDR[['time']] <- c('biSBM'=bisbm.timing,'SBM'=sbm.timing)

print(FDRTDR)


saveRDS(FDRTDR, file=paste0(OutFolder,'fdrtdr_',file.end,'_',jid,'.rds'))

