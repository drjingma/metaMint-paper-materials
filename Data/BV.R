## script to perform stability analysis on the BV data set. 

rm(list=ls())
source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")
source("biSBM/stability_selection.R")

## Toggle between lines 12 and 13 to run with model selection or for a specific model
# args <- c('Gauss0Var1','0','pearson','4','1') # run the method on a specific model with Q1=4 and Q2 = 1
args = commandArgs(trailingOnly = T) # run the method with parameters specified in the bash file
print(args)

output_dir <- "../Data/BV_Analysis/"
d <- readRDS(paste0(output_dir, 'BV_teststatistics_correlation.rds'))
model <- as.character(args[1])
exclusive_Rows_Cols <- as.integer(args[2])
type <- as.character(args[3])
Q1 <- as.integer(args[4])
Q2 <- as.integer(args[5])
ICL <- (Q1*Q2>0)
Q1_range <- switch( ICL + 1,  # because switch() in numeric mode is 1-indexed
                    seq_len(10),      # case ICL == 0
                    Q1                # case ICL == 1
)
  
Q2_range <-  switch( ICL + 1,  # because switch() in numeric mode is 1-indexed
                     seq_len(10),      # case ICL == 0
                     Q2                # case ICL == 1
)

file.path <- paste0(output_dir, model)
if (!dir.exists(file.path)){
  dir.create(file.path, showWarnings = T, recursive = TRUE)
  print(file.path)
}
file.ending <- paste0('BV_differential_exclusive', exclusive_Rows_Cols, '_',type, '_Q1',Q1,'_Q2',Q2)

dataMatrix <- switch(type,
                     perm = {d$perm},
                     cl = {d$cl},
                     asym = {d$asym},
                     pearson = {d$pearson})

stab <- list()
for (Q1 in Q1_range){
  for (Q2 in Q2_range){
    key <- sprintf("Q1=%d_Q2=%d", Q1, Q2)
    cat("Running Q1 =", Q1, "Q2 =", Q2, "\n")
    stab[[key]] <- stability_nobiSBM(
      dataMatrix,
      model = 'Gauss01',
      sbmSize = list(Q1 = Q1, Q2 = Q2),
      exclusive_Rows_Cols = FALSE,
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
    
    saveRDS(stab, file = paste0(file.path,'/', file.ending,'_Jaccard_stability.rds'))
    
  }
}


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

xtable::xtable(x = scores_df %>% arrange(desc(norm_stab))  %>% head(10))
