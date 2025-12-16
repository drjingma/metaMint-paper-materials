rm(list=ls())

source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")
source("biSBM/stability_selection.R")
library(dplyr)
library(xtable)

args <- c("1",'100', "150","200","3","3","1","1","1","3","0.1","0.8",'Gauss','Gauss01', '1')
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

file.ending <- paste0("/stab_n1_",n1,"_n2_",n2,"_Q1_",Q1,"_Q2_",Q2, '_alt', mu_diag, mu_offdiag, '_sd', nu0[2], nu, '_rep', run_rep, '_nreps_', nreps, "_ICL",ICL, ".rds")

if (ICL){
  Q1_range <- seq_len(5)
  Q2_range <- seq_len(5)
} else {
  Q1_range = Q1
  Q2_range = Q2
}
stab <- list()
for (Q1 in Q1_range){
  for (Q2 in Q2_range){
    key <- sprintf("Q1=%d_Q2=%d", Q1, Q2)
    cat("Running Q1 =", Q1, "Q2 =", Q2, "\n")
    stab[[key]] <- # Then call
      stability_nobiSBM(
        d$dataMatrix,
        model = 'Gauss01',
        sbmSize = list(Q1 = Q1, Q2 = Q2),
        exclusive_Rows_Cols = FALSE,
        num.iter = 5,
        initParam = list(nbOfbeta = NULL, nbOfPointsPerbeta = NULL,
                         maxNbOfPasses = NULL, minNbOfPasses = 1),
        B = 50,
        row_frac = 0.8,
        col_frac = 0.8,
        select_by = "first",
        fix_Q1 = NULL,
        fix_Q2 = NULL,
        seed = NULL,
        keep_per_run = FALSE,
        verbose = TRUE
      )
  }
}
saveRDS(stab, file = paste0(output_folder, file.ending))

file.ending <- "/stab_n1_150_n2_200_Q1_3_Q2_3_alt13_sd11_rep1_nreps_100_ICL1.rds"
stab <- readRDS(paste0(output_folder, file.ending))
scores <- t(sapply(stab, function(s) {
  bic_df <- get_bicluster_consensus_from_stab(s)
  
  # random baseline + normalized bicluster stability
  bic_df$c_rand <- compute_random_bic_cons(bic_df, s)
  bic_df$c_norm <- (bic_df$bic_cons - bic_df$c_rand) / (1 - bic_df$c_rand)
  
  # model-level normalized stability (size-weighted)
  model_norm_stab <- with(
    bic_df,
    weighted.mean(c_norm, w = n_rows * n_cols, na.rm = TRUE)
  )
  
  # grab ICL from the full-fit object
  model_ICL <- s$full_fit$sbmParam$ICL  # adjust if your ICL is stored differently
  
  # return both as a named vector
  c(norm_stab = model_norm_stab,
    ICL       = model_ICL)
}))

scores_df <- as.data.frame(scores)

# scores_df <- data.frame(S_weighted = scores)
scores_df %>% arrange(desc(norm_stab))  %>% head(10)

scores_df %>% arrange(desc(ICL))  %>% head(10)

xtable(x = scores_df %>% arrange(desc(norm_stab))  %>% head(10))
