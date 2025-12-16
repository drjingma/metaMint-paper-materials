rm(list=ls())

source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")
source("biSBM/FDR-R-Code/adaptZ.func.R")
library(dplyr)

args = commandArgs(trailingOnly = T)
print(args)

run_rep <- as.integer(args[1])
nreps <- as.integer(args[2])
n1 <- as.integer(args[3])
n2 <- as.integer(args[4])
Q1 <- as.integer(args[5])
Q2 <- as.integer(args[6])
nu0 <- c(0,as.numeric(args[7]))
nu <- as.numeric(args[8])
mu_diag <- as.numeric(args[9])
mu_offdiag <- as.numeric(args[10])
prob <- c(as.numeric(args[11]),
          as.numeric(args[12]))
modelFamily <- as.character(args[13])
model <- as.character(args[14])
ICL <- as.integer(args[15]) # ICL==2: indicating misleading the model to incorrect B1,B2


exclusive_Rows_Cols <- F

output_folder = paste0("../Output/newResults_biSBM/",modelFamily,'/', model, "/", paste0(prob, collapse = "_"))

if (!dir.exists(output_folder)){
  dir.create(output_folder, showWarnings = T, recursive = TRUE)
  print(output_folder)
}

file.ending <- paste0("/data_n1_",n1,"_n2_",n2,"_Q1_",Q1,"_Q2_",Q2, '_alt', mu_diag, mu_offdiag, '_sd', nu0[2], nu, '_rep', run_rep, '_nreps_', nreps,".rds")
check_file = paste0(output_folder, file.ending)

if_have_file = tryCatch(
  {
    d <- readRDS(check_file)
  }
  ,error =  function(e) e
)

if ("simpleError" %in% class(if_have_file)){
  theta <- list(alpha1 = rep(1/Q1,Q1),
                alpha2 = rep(1/Q2,Q2),
                nu0 = c(0,1))
  theta$nu <- list(mean = matrix(mu_offdiag,Q1,Q2)  %>% `diag<-`(., mu_diag),
                   sd = matrix(1,Q1,Q2))
  theta$pi <- matrix(prob[1],Q1,Q2) %>% `diag<-`(., prob[2])
  
  d <- rnbisbm(n1,n2,theta = theta, modelFamily = modelFamily)
  saveRDS(d, file= paste0(output_folder, file.ending))
}

file.ending <- paste0("/res_n1_",n1,"_n2_",n2,"_Q1_",Q1,"_Q2_",Q2, '_alt', mu_diag, mu_offdiag, '_sd', nu0[2], nu, '_rep', run_rep, '_nreps_', nreps, "_ICL",ICL, ".rds")

A <- d$latentAdj
dataVector <- as.vector(d$dataMatrix)
pvals <- 2*pnorm(abs(dataVector),lower.tail = FALSE)

start.time = Sys.time()
res <- fitNobiSBM(d$dataMatrix, model=model,
                  exclusive_Rows_Cols = F, 
                  sbmSize = switch(
                    ICL + 1,  # because switch() in numeric mode is 1-indexed
                    list(Q1 = Q1, Q2 = Q2, explor = 1.5),        # case ICL == 0
                    list(Q1 = 1:5, Q2 = 1:5, explor = 1.5),      # case ICL == 1
                    list(Q1 = 1, Q2 = 1, explor = 1.5),          # case ICL == 2
                    list(Q1 = 5, Q2 = 5, explor = 1.5)           # case ICL == 3
                  ),
                  initParam = list(nbOfbeta=1, nbOfPointsPerbeta=NULL,
                                   maxNbOfPasses=2, minNbOfPasses=1),
                  nbCores = 1)
end.time = Sys.time()
bisbm.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

best.id <- which(sapply(res,function(a) a$sbmParam$ICL)==getBestQ(res, exclusive_Rows_Cols)$ICL)
currentSolution <- res[[best.id]]
currentSolution$theta
alpha.range <- c(0.005,0.025,0.05,0.1,0.15,0.25)
gnbiSBM <- graphInferenceNobiSBM(d$dataMatrix, 
                                 currentSolution$clustering_row, 
                                 currentSolution$clustering_col, 
                                 currentSolution$theta,alpha = 0.05)
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

saveRDS(FDRTDR, file= paste0(output_folder, file.ending))



