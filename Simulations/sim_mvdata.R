# simulate multi-view data from a given covariance matrix and mean
# Two-sample case

rm(list=ls())

library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(pheatmap)
library(phyloseq)

source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")
source("biSBM/FDR-R-Code/adaptZ.func.R")

# args ----
# args <- c("1",'100', "49","50","3","3","400","4","1","0","0.1","0.8",'ln','Gauss01', '0')
args = commandArgs(trailingOnly = T)
print(args)

run_rep <- as.integer(args[1])
nreps <- as.integer(args[2])
n1 <- as.integer(args[3])
n2 <- as.integer(args[4])
Q1 <- as.integer(args[5])
Q2 <- as.integer(args[6])
m2 <- m1 <- as.numeric(args[7])
pseudocount <- as.numeric(args[8])
use_clr <- as.numeric(args[9])
even_depth <- as.numeric(args[10])
prob <- c(as.numeric(args[11]),
          as.numeric(args[12]))
modelFamily <- as.character(args[13])
model <- as.character(args[14])
ICL <- as.integer(args[15])


shrink_factor <- 0.1
show_plot <- F
exclusive_Rows_Cols <- F
alpha.range <- c(0.005,0.025,0.05,0.1,0.15,0.25)

output_folder = paste0("../Output/newResults_biSBM/",modelFamily,'/', model, "/mvdata")
if (!dir.exists(output_folder)){
  dir.create(output_folder, showWarnings = T, recursive = TRUE)
  print(output_folder)
}

file.ending <- paste0("/data_n1_",n1,"_n2_",n2, "_Q1_",Q1, "_Q2_",Q2, '_m', m1,  
                      "_clr", use_clr, 
                      '_c', pseudocount,
                      '_prob', paste0(prob, collapse = "_"), '_rep', run_rep, '_nreps_', nreps,".rds")

d <- readRDS(paste0("../Output/newResults_biSBM/negbin/", model,
                         "/correlation_n1_",n1,"_n2_",n2, "_Q1_",Q1, "_Q2_",Q2, '.rds'))
Sigma1 <- d$Sigma1
Sigma2 <- d$Sigma2
load("../Data/BV_Analysis/BV.rda")
ASV <- as(otu_table(BV1), "matrix")
metabolite <- t(dataset$dat)

if (modelFamily == 'negbin') {
dat1 <- synth_comm_from_XY(X = t(ASV+pseudocount)[,1:n1], 
                           Y = exp(metabolite[,1:n2]),
                           n = m1,
                           Sigma = Sigma1, 
                           library_size = NULL,
                           return_latent = TRUE,
                           library_scale = switch(
                             even_depth + 1,  # because switch() in numeric mode is 1-indexed
                             round(runif(m1,min=1, max=10)),        # case even_depth == 0
                             rep(1,m1)                              # case even_depth == 1
                           ),
                           seed = NULL)

dat2 <- synth_comm_from_XY(X = t(ASV+pseudocount)[,1:n1], 
                           Y = exp(metabolite[,1:n2]),
                           n = m2,
                           Sigma = Sigma2, 
                           return_latent = TRUE,
                           library_scale = switch(
                             even_depth + 1,  # because switch() in numeric mode is 1-indexed
                             round(runif(m2,min=1, max=10)),        # case even_depth == 0
                             rep(1,m2)                              # case even_depth == 1
                           ),
                           seed = NULL)
}

if (modelFamily == 'lnm'){
  mu = runif(n1+n2, 0, 4)
  D = diag(runif(n1+n2,1,3))
  Sigma1 <- D %*% Sigma1 %*% D 
  Sigma2 <- D %*% Sigma2 %*% D 
  dat1 <- para_data_generation(m1, n1, mu, Sigma1, library_scale = switch(
    even_depth + 1,  # because switch() in numeric mode is 1-indexed
    sample(rowSums(ASV), m1, replace = TRUE),        # case even_depth == 0
    rep(1e6,m1)                                      # case even_depth == 1
  ))
  dat2 <- para_data_generation(m2, n1, mu, Sigma2, library_scale = switch(
    even_depth + 1,  # because switch() in numeric mode is 1-indexed
    sample(rowSums(ASV), m2, replace = TRUE),        # case even_depth == 0
    rep(1e6,m2)                                      # case even_depth == 1
  ))
}

if (modelFamily == 'ln'){
  mu = runif(n1+n2, -1*pseudocount, 4)
  D = diag(runif(n1+n2,1,3))
  Sigma1 <- D %*% Sigma1 %*% D 
  Sigma2 <- D %*% Sigma2 %*% D 
  dat1 <- para_data_generation_ln(m1, n1, mu, Sigma1)
  dat2 <- para_data_generation_ln(m2, n1, mu, Sigma2)
}
print(mean(dat1$X==0))

sample1 <- list(X = t(dat1$Z[,1:n1]), Y = t(dat1$Z[,-(1:n1)]))
sample2 <- list(X = t(dat2$Z[,1:n1]), Y = t(dat2$Z[,-(1:n1)]))

if (use_clr){
  sample1$X <- t(malabutils::clr.epsilon(dat1$X, const = 0.01))
  sample2$X <- t(malabutils::clr.epsilon(dat2$X, const = 0.01))
  sample1$Y <- t(log(dat1$Y))
  sample2$Y <- t(log(dat2$Y))
}
# if (use_fisher){
#   R <- equal_cor_matrix_test(dat1$Z, dat2$Z, alternative = 'two.sided', p_adjust = 'none')
#   dataMatrix <- R$z_matrix[1:n1, -(1:n1)]
# } else {
dataMatrix <- teststat_general_2sample(sample1,sample2)
# }
saveRDS(dataMatrix, 
        file= paste0(output_folder, file.ending))

# }

A <- d$A
dataVector <- as.vector(dataMatrix)
pvals <- 2*pnorm(abs(dataVector),lower.tail = FALSE)

if (show_plot){
  print(range(dataMatrix))
  pheatmap(dataMatrix[order(d$latentZ[[1]]), order(d$latentZ[[2]])],
           cluster_rows = F, cluster_cols = F)
  # Plot histogram, scaled to density
  hist(dataVector, breaks = 30, freq = FALSE,
       col = "gray80", border = "white",
       main = "Histogram of z with Standard Normal Overlay",
       xlab = "z-values")
  # Overlay standard normal density
  curve(dnorm(x), col = "red", lwd = 2, add = TRUE)
}

res <- fitNobiSBM(dataMatrix, model=model, 
                  exclusive_Rows_Cols = exclusive_Rows_Cols,
                  sbmSize = switch(
                    ICL + 1,  # because switch() in numeric mode is 1-indexed
                    list(Q1 = Q1, Q2 = Q2, explor = 1.5),        # case ICL == 0
                    list(Q1 = 1:5, Q2 = 1:5, explor = 1.5)       # case ICL == 1
                  ),
                  initParam = list(nbOfbeta=1, nbOfPointsPerbeta=NULL,
                                   maxNbOfPasses=2, minNbOfPasses=1),
                  nbCores = 1)
best.id <- which(sapply(res,function(a) a$sbmParam$ICL)==getBestQ(res, exclusive_Rows_Cols)$ICL)
currentSolution <- res[[best.id]]
gnbiSBM <- graphInferenceNobiSBM(dataMatrix,
                                 currentSolution$clustering_row,
                                 currentSolution$clustering_col,
                                 currentSolution$theta, alpha=0.05)

prep <- adaptZ_prepare(dataVector, gamma = 0.1, model = "Gauss01")
get_metric <- function(a) {
  mat <- matrix(0, 2, 4)
  rownames(mat) <- c("FDR",'TDR')
  colnames(mat) <- c('New','BH','Storey','SC')
  Ahat_bisbm <- matrix(1*(gnbiSBM$qvalues<a), n1, n2)
  FDR <- sum(Ahat_bisbm[A==0])/ max(sum(Ahat_bisbm),1)
  TDR <- sum(Ahat_bisbm[A==1])/ sum(A==1)
  mat[,1] <- c(FDR=FDR,TDR=TDR)
  
  # BH
  qval_results <- p.adjust(pvals,"BH")
  Ahat_bh <- matrix(1*(qval_results<a), n1, n2)
  FDR <- sum(Ahat_bh[A==0])/ max(sum(Ahat_bh),1)
  TDR <- sum(Ahat_bh[A==1])/ sum(A==1)
  mat[,2] <- c(FDR=FDR,TDR=TDR)
  
  # Storey's q
  qval_results <- qvalue::qvalue(pvals,fdr.level=NULL)$qvalues
  Ahat_q <- matrix(1*(qval_results<a), n1, n2)
  FDR <- sum(Ahat_q[A==0])/ max(sum(Ahat_q),1)
  TDR <- sum(Ahat_q[A==1])/ sum(A==1)
  mat[,3] <- c(FDR=FDR,TDR=TDR)
  
  # SC07
  adaptiveZ <- adaptZ_apply(prep,a)
  
  Ahat_sc <- matrix(0,n1,n2)
  Ahat_sc[adaptiveZ$re] <- 1
  FDR <- sum(Ahat_sc[A==0])/ max(sum(Ahat_sc),1)
  TDR <- sum(Ahat_sc[A==1])/ sum(A==1)
  mat[,4] <- c(FDR=FDR,TDR=TDR)
  
  # Compare overlap with the new method
  # fraction of discoveries in each benchmark
  prop <- sapply(list(Ahat_bh, Ahat_q, Ahat_sc), function(m) {
    overlap <- length(intersect(which(Ahat_bisbm==1), which(m==1)))
    overlap / max(length(which(m==1)),1)   
  })     
  
  
  list(fdr = mat,
       prop = prop)
}

names(alpha.range) <- as.character(alpha.range)
FDRTDR <- lapply(alpha.range, get_metric)
FDRTDR$Q <- getBestQ(res,F)$Q1
FDRTDR$randIndex <- c('row' = fossil::adj.rand.index(currentSolution$clustering_row, d$latentZ[[1]]),
                      'col' = fossil::adj.rand.index(currentSolution$clustering_col, d$latentZ[[2]]))
# print(FDRTDR)
file.ending <- paste0("/res_n1_",n1,"_n2_",n2,"_Q1_",Q1,"_Q2_",Q2, '_m', m1, 
                      '_clr', use_clr, 
                      '_c', pseudocount,  
                      '_even', even_depth,
                      '_prob', paste0(prob, collapse = "_"), 
                      '_rep', run_rep, '_nreps_', nreps, "_ICL", ICL, ".rds")

saveRDS(FDRTDR, file= paste0(output_folder, file.ending))


