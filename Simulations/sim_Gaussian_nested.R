
source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")
source("biSBM/FDR-R-Code/adaptZ.func.R")
library(dplyr)

jid=commandArgs(trailingOnly=T)[1]
jid=as.numeric(jid)

OutFolder <- "../Output/biSBM/"

today <- "20250517"
model <- "Gauss01"
setting <- 4

n1 <- 150
n2 <- 200

Q1 <- 2
Q2 <- 2
mu <- 2

exclusive_Rows_Cols <- F#ifelse(Q1==Q2, T, F)

file.end <- paste0(today,"_n1_",n1,"_n2_",n2,"_Q1_",Q1,"_Q2_",Q2, '_',model, '_setting', setting, '_mu')

A <- matrix(0, n1, n2)
A[1,] <- 1
A[,1] <- 1
dataMatrix <- matrix(rnorm(n1*n2), n1, n2)
dataMatrix[A==1] <- rnorm(sum(A==1), mean = mu)

A.full <- rbind(cbind(matrix(0,n1,n1),A),
                cbind(t(A), matrix(0,n2,n2)))

R <- rbind(cbind(matrix(rnorm(n1^2),n1,n1),dataMatrix),
           cbind(t(dataMatrix), matrix(rnorm(n2^2),n2,n2)))

dataVector <- as.vector(dataMatrix)
pvals <- 2*pnorm(abs(dataVector),lower.tail = FALSE)

## fit biSBM ----
start.time = Sys.time()
res <- fitNobiSBM(dataMatrix, model=model,
                  exclusive_Rows_Cols = F, 
                  # sbmSize = list(Q1=Q1, Q2=Q2, explor=1.5),
                  sbmSize = list(Q1=1:5, Q2=1:5, explor=1.5),
                  initParam = list(nbOfbeta=1, nbOfPointsPerbeta=NULL,
                                   maxNbOfPasses=2, minNbOfPasses=1),
                  nbCores = 1)
end.time = Sys.time()
bisbm.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

best.id <- which(sapply(res,function(a) a$sbmParam$ICL)==getBestQ(res, exclusive_Rows_Cols)$ICL)
currentSolution <- res[[best.id]]
alpha.range <- c(0.005,0.025,0.05,0.1,0.15,0.25)
FDRTDR <- list()
gnbiSBM <- graphInferenceNobiSBM(dataMatrix, 
                                 currentSolution$clustering_row, 
                                 currentSolution$clustering_col, 
                                 currentSolution$theta,alpha = 0.05)
FDRTDR[['biSBM']] <- sapply(alpha.range, function(a) {
  Ahat <- matrix(1*(gnbiSBM$qvalues<a), n1, n2)
  FDR <- sum(Ahat[A==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A==1])/ sum(A==1)
  data.frame(FDR=FDR,TDR=TDR)
})

## noisySBM ----
library(noisySBM)
start.time = Sys.time()
res.sbm <- fitNSBM(R, model = model,
                   sbmSize = list(Qmin = 1, Qmax = Q1+Q2+1, explor = 1.5),
                   # sbmSize = list(Qmin = Q1+Q2, Qmax = Q1+Q2, explor = 1.5),
                   initParam = list(nbOfTau = 1, nbOfPointsPerTau = NULL,
                                    maxNbOfPasses = 2,
                                    minNbOfPasses = 1))
end.time = Sys.time()
sbm.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

best.id <- which(sapply(res.sbm,function(a) a$sbmParam$ICL)==noisySBM::getBestQ(res.sbm)$ICL)
currentSolution.sbm <- res.sbm[[best.id]]

saveRDS(R, file= paste0(OutFolder, file.end,'_data_',jid,'.rds'))
saveRDS(res.sbm, file= paste0(OutFolder, file.end,'_SBMresults_',jid,'.rds'))

gnSBM <- noisySBM::graphInference(R,
                                  currentSolution.sbm$clustering,
                                  currentSolution.sbm$theta, alpha=0.05)

FDRTDR[['SBM']] <- sapply(alpha.range, function(a) {
  Ahat <- matrix(0, n1+n2, n1+n2)
  Ahat[lower.tri(Ahat)] <- (gnSBM$qvalues < a)
  Ahat <- Ahat + t(Ahat)
  FDR <- sum(Ahat[A.full==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A.full==1])/ sum(A.full==1)
  data.frame(FDR=FDR,TDR=TDR)
})

FDRTDR[['SBM.2']] <- sapply(alpha.range, function(a) {
  Ahat <- matrix(0, n1+n2, n1+n2)
  Ahat[lower.tri(Ahat)] <- (gnSBM$qvalues < a)
  Ahat <- Ahat + t(Ahat)
  Ahat <- Ahat[1:n1, -(1:n1)]
  FDR <- sum(Ahat[A==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A==1])/ sum(A==1)
  data.frame(FDR=FDR,TDR=TDR)
})

FDRTDR[['BH']] <- sapply(alpha.range, function(a) {
  qval_results <- p.adjust(pvals,"BH")
  Ahat <- matrix(1*(qval_results<a), n1, n2)
  FDR <- sum(Ahat[A==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A==1])/ sum(A==1)
  data.frame(FDR=FDR,TDR=TDR)
})
FDRTDR[['Storey']] <- sapply(alpha.range, function(a) {
  qval_results <- qvalue::qvalue(pvals,fdr.level=NULL)$qvalues
  Ahat <- matrix(1*(qval_results<a), n1, n2)
  FDR <- sum(Ahat[A==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A==1])/ sum(A==1)
  data.frame(FDR=FDR,TDR=TDR)
})

FDRTDR[['SC07']] <- sapply(alpha.range, function(a) {
  adaptiveZ <- adaptZ.func(dataVector, a, 0.1, 'Gauss01')
  # the threshold
  threshold <- adaptiveZ$th
  # number of rejected hypotheses
  k <- adaptiveZ$nr
  # the rejected hypotheses
  rh <- adaptiveZ$re
  
  Ahat <- matrix(0,n1,n2)
  Ahat[adaptiveZ$re] <- 1
  FDR <- sum(Ahat[A==0])/ max(sum(Ahat),1)
  TDR <- sum(Ahat[A==1])/ sum(A==1)
  data.frame(FDR=FDR,TDR=TDR)
})
FDRTDR$Q <- getBestQ(res)$Q1

FDRTDR$time <- c('bisbm' = bisbm.timing, 'sbm' = sbm.timing)

saveRDS(FDRTDR, file= paste0(OutFolder, file.end,'_',jid,'.rds'))


