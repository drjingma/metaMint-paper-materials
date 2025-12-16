source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")

args <- c("1",'100', "49","50","3","3","250","0","1","0","0.1","0.8",'negbin','Gauss01', '0')

n1 <- as.integer(args[3])
n2 <- as.integer(args[4])
Q1 <- as.integer(args[5])
Q2 <- as.integer(args[6])
modelFamily <- as.character(args[13])
model <- as.character(args[14])

load("../Data/BV_Analysis/BV.rda")
ASV <- as(otu_table(BV1), "matrix")
# ASV <- ASV[sample(1:25, n1, replace = T),]
ASV_clr <- malabutils::clr.epsilon(t(ASV), const = 0.01)
metabolite <- t(dataset$dat)
theta <- list(alpha1 = rep(1/Q1,Q1),
              alpha2 = rep(1/Q2,Q2), 
              nu0 = c(0,1))
theta$nu <- list(mean = matrix(3,Q1,Q2),
                 sd = matrix(1,Q1,Q2)) 
theta$pi <- matrix(prob[1],Q1,Q2) %>% `diag<-`(., prob[2])

# generate the probability matrix
d <- rnbisbm(n1, n2, theta = theta)

# per-(block,block) multipliers
rho_by_block <- matrix(c(1,-1,1,
                         -1,-1,-1,
                         1,-1,1), 
                       nrow = max(d$latentZ[[1]]), byrow = TRUE)
out <- build_corr_from_bipartite(
  d$latentAdj, d$latentZ[[1]], 
  Zy=d$latentZ[[2]], 
  ASV_clr[,1:n1], 
  metabolite[,1:n2],
  cross_rho = 1,
  rho_by_block = rho_by_block,
  method_within = "shrink", shrink_lambda = 0.5,
  corr_type = "pearson",
  jitter_xy = 0.0,
  project_nearPD = TRUE,
  target_spec = 0.99,
  cross_mode = 'whitened'
)

Sigma1 <- out$R
Sigma2 <- Sigma1
Sigma2[1:n1, -(1:n1)] <- - Sigma1[1:n1, -(1:n1)] 
Sigma2[-(1:n1), (1:n1)] <- - Sigma1[-(1:n1), (1:n1)]

saveRDS(list(Sigma1 = Sigma1, Sigma2 = Sigma2, A = d$latentAdj, latentZ = d$latentZ),
        file = paste0("../Output/newResults_biSBM/",modelFamily,'/', model,
                      "/correlation_n1_",n1,"_n2_",n2, "_Q1_",Q1, "_Q2_",Q2, '.rds'))
