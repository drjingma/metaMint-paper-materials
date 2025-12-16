#' transform a pair of nodes (i,j) into an identifying integer
#'
#' Associates an identifying integer with a pair of nodes (i,j)
#'
#' @details returns the row number of the matrix build by listNodePairs(n)
#'     containing the pair (i,j)
#'
#' @param i scalar or vector
#' @param j scalar or vector, same length as i
#' @param n number of vertices
#' @param directed booelan to indicate whether the model is directed or undirected
# convertNodePair <- function(i,j,n, directed){
#   if (sum((i>n) | (j>n))>0){
#     stop("Your index is out of range")
#   }
#   if (directed){#directed case
#     dyads <- (i-1)*(n-1)+j-(i<j)
#   } else {#undirected case
#     dyads <- c(0,cumsum((n-1):1))[pmin(i,j)] + abs(j-i)
#   }
#   return(dyads)
# }



#' returns a list of all possible node pairs (i,j)
#'
#' @param n1 number of nodes on side one
#' @param n2 number of nodes on side two
#' @param directed indicates if the graph is directed
#'
#' @return a 2-column matrix with all possible node pairs (i,j)
listNodePairs <- function(n1, n2, directed=FALSE){
  return(expand.grid(1:n1,1:n2))
}


#' transform a pair of block identifiers (q,l) into an identifying integer
#'
#' this is the inverse function of convertGroupPairIdentifierBipartite()
#'
#' @param q indicator of a latent block
#' @param l indicator of a latent block; can be a vector
#' @param Q1 number of latent blocks on side one
#' @param Q2 number of latent blocks on side two
#' @param directed indicates if the graph is directed
convertGroupPair <- function(q, l, Q1, Q2){
  all <- listNodePairs(Q1,Q2)
  index <- match(sapply(l, function(a) paste(q, a, collapse =" ")),apply(all, 1, paste, collapse =" "))
  return(index)
}


#' takes a scalar index of a group pair (q,l) and returns the values q and l
#'
#' this is the inverse function of convertGroupPair()
#'
#' @param ind_ql indicator for a pair of latent blocks
#' @param Q1 number of latent blocks on side one
#' @param Q2 number of latent blocks on side two
# convertGroupPairIdentifierBipartite <- function(ind_ql, Q1, Q2){
#   q <- ind_ql %% Q1
#   if (q==0){
#     q <- Q1
#   }
#   l <- ceiling(ind_ql/Q1)
#   return(c(q,l))
# }

#' takes a scalar indice of a group pair (q,l) and returns the values q and l
#'
#' this is the inverse function of convertGroupPair()
#'
#' @param ind_ql indicator for a pair of latent blocks
#' @param Q number of latent blocks
convertGroupPairIdentifier <- function(ind_ql, Q){
  w <- cumsum((Q-1):1)
  q <- which.max(ind_ql<=w)
  w <- c(0, w)
  l <- ind_ql - w[q] + q
  return(c(q,l))
}


#' corrects values of the variational parameters tau that are too close to the 0 or 1
#'
#' @param tau variational parameters
correctBeta <- function(tau){
  tau <- pmin(tau,.Machine$double.xmax)
  tau <- pmax(tau,.Machine$double.xmin)
  tau <- tau/sum(tau)
  tau <- pmin(tau,1-1e-7)
  tau <- pmax(tau,1e-7)
  tau <- tau/sum(tau)
  
  return(tau)
}


#' evaluate the density in the current model
#'
#' @param x vector with points where to evaluate the density
#' @param nu distribution parameter
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}, \code{t}
modelDensity <- function(x, nu, modelFamily='Gauss'){
  if (modelFamily=='Gauss')
    res <- stats::dnorm(x, nu[1], nu[2])
  if (modelFamily=='Gamma')
    res <- stats::dgamma(x, nu[1], nu[2])
  if (modelFamily=='t')
    res <- stats::dt(x, ncp = nu[1], df = nu[2])
  res[res<=.Machine$double.eps] <- .Machine$double.eps
  return(res)
}


#' Evaluate beta_q_1*beta_l_2 in the noisy bipartite stochastic block model
#'
#' @param q indicator of a latent block on side one
#' @param l indicator of a latent block on side two
#' @param beta variational parameters in a list of length 2
getBetaql <- function(q, l, beta){
  beta1 <- beta[[1]]
  beta2 <- beta[[2]]
  n1 <- ncol(beta1)
  n2 <- ncol(beta2)
  Q1 <- nrow(beta1)
  Q2 <- nrow(beta2)
  
  # would like to get a matrix output
  if (Q1 == 1 && Q2 == 1) # one block on both sides
    betaql <- matrix(1,n1,n2)
  else{
    betaql <- outer(beta1[q, ],beta2[l, ], '*')
  }
  return(betaql)
}

#' Function to calculate test statistics
#' @param X,Y feature by sample matrices
teststat_general <- function(X,Y=NULL){
  n1 <- ncol(X)
  if (is.null(Y)){
    Y <- X
  }
  empcor <- cor(t(X),t(Y),method = 'pearson')
  X.scale <- t(scale(t(X), center = T, scale = T))
  Y.scale <- t(scale(t(Y), center = T, scale = T))
  Theta <- matrix(0, nrow(X), nrow(Y))
  for (i in 1:nrow(Theta)){
    for (j in 1:ncol(Theta)){
      Theta[i,j] <- mean((2 * X.scale[i,] * Y.scale[j,] - empcor[i,j] * X.scale[i,]^2 - empcor[i,j] * Y.scale[j,]^2)^2)
    }
  }
  dataMatrix <- 2*empcor*sqrt(n1/Theta)
  
  return(list(mean = empcor, sd = Theta, statistic = dataMatrix))
}
#' Function to calculate two-sample test statistics 
teststat_general_2sample <- function(sample1, sample2){
  n1 <- ncol(sample1$X)
  n2 <- ncol(sample2$X)
  
  out1 <- teststat_general(sample1$X,sample1$Y)
  out2 <- teststat_general(sample2$X,sample2$Y)
  
  dataMatrix <- 2*(out1$mean - out2$mean)/sqrt(out1$sd/n1 + out2$sd/n2)
  
  dataMatrix
}

# Simultaneous two-sample test of equal correlations (MVN, independent groups)
# H0: R_X[i,j] = R_Y[i,j] for all i<j
equal_cor_matrix_test <- function(X, Y,
                                  alternative = c("two.sided","less","greater"),
                                  p_adjust = c("BH","bonferroni","holm","BY","none")) {
  alternative <- match.arg(alternative)
  p_adjust    <- match.arg(p_adjust)
  
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (ncol(X) != ncol(Y)) stop("X and Y must have the same number of columns (variables).")
  p  <- ncol(X)
  
  # containers
  r1 <- r2 <- matrix(NA_real_, p, p)
  n1 <- n2 <- matrix(NA_integer_, p, p)
  z  <- pval <- matrix(NA_real_, p, p)
  
  # pairwise complete cases per pair (allows missingness)
  for (i in 1:(p-1)) for (j in (i+1):p) {
    cc1 <- stats::complete.cases(X[, i], X[, j])
    cc2 <- stats::complete.cases(Y[, i], Y[, j])
    n1_ij <- sum(cc1); n2_ij <- sum(cc2)
    n1[i,j] <- n1_ij; n2[i,j] <- n2_ij
    if (n1_ij > 3 && n2_ij > 3) {
      r1[i,j] <- stats::cor(X[cc1, i], X[cc1, j])
      r2[i,j] <- stats::cor(Y[cc2, i], Y[cc2, j])
      z1 <- atanh(max(min(r1[i,j], 0.999999), -0.999999))
      z2 <- atanh(max(min(r2[i,j], 0.999999), -0.999999))
      se <- sqrt(1/(n1_ij - 3) + 1/(n2_ij - 3))
      z[i,j] <- (z1 - z2) / se
      # p-values
      pval[i,j] <- switch(alternative,
                          two.sided = 2 * pnorm(-abs(z[i,j])),
                          greater   = 1 - pnorm(z[i,j]),   # H1: r1 > r2
                          less      = pnorm(z[i,j]))       # H1: r1 < r2
    }
  }
  
  # mirror to full matrices for convenience
  r1 <- r1 + t(r1); diag(r1) <- 1
  r2 <- r2 + t(r2); diag(r2) <- 1
  # Make symmetric matrices without clobbering finite values
  z[lower.tri(z)]    <- t(z)[lower.tri(z)]
  pval[lower.tri(pval)] <- t(pval)[lower.tri(pval)]
  
  # vectorized upper-tri results
  ut <- which(upper.tri(z), arr.ind = TRUE)
  z_ut    <- z[ut]
  p_ut    <- pval[ut]
  n1_ut   <- n1[ut]; n2_ut <- n2[ut]
  r1_ut   <- r1[ut]; r2_ut <- r2[ut]
  
  # multiple testing adjustment (on the finite p-values)
  ok <- is.finite(p_ut)
  padj <- rep(NA_real_, length(p_ut))
  if (any(ok)) padj[ok] <- p.adjust(p_ut[ok], method = p_adjust)
  
  # Global tests
  m <- sum(is.finite(z_ut))                   # number of tested edges
  z2_sum <- sum(z_ut[is.finite(z_ut)]^2)
  p_global_chisq <- if (m > 0) pchisq(z2_sum, df = m, lower.tail = FALSE) else NA_real_
  
  zmax <- max(abs(z_ut[is.finite(z_ut)]), na.rm = TRUE)
  p_global_maxbonf <- if (is.finite(zmax)) min(1, (2 * (1 - pnorm(zmax))) * m) else NA_real_
  
  # Assemble a tidy table
  res_table <- data.frame(
    i = ut[,1], j = ut[,2],
    n1 = n1_ut, n2 = n2_ut,
    r1 = r1_ut, r2 = r2_ut,
    z  = z_ut,
    p  = p_ut,
    padj = padj,
    stringsAsFactors = FALSE
  )
  
  structure(list(
    method = "Simultaneous two-sample Fisher z-tests for equality of correlations (MVN)",
    alternative = alternative,
    p_adjust = p_adjust,
    table = res_table[order(res_table$padj, na.last = TRUE), ],
    z_matrix = z,
    p_matrix = pval,
    r1 = r1, r2 = r2,
    global = list(
      m_edges = m,
      chisq_sum_z2 = z2_sum,
      df = m,
      p_global_chisq = p_global_chisq,
      zmax = zmax,
      p_global_max_bonf = p_global_maxbonf
    )
  ), class = "equal_cor_matrix_test")
}

# Pretty print
print.equal_cor_matrix_test <- function(x, ...) {
  cat(x$method, "\n")
  cat("Alternative:", x$alternative, "| p-adjust:", x$p_adjust, "\n")
  g <- x$global
  cat(sprintf("Global χ² (sum z^2): %.2f on %d df, p=%.3g\n",
              g$chisq_sum_z2, g$df, g$p_global_chisq))
  cat(sprintf("Global max |z|: %.2f, Bonferroni p=%.3g\n",
              g$zmax, g$p_global_max_bonf))
  cat("\nTop edges (by adjusted p):\n")
  print(utils::head(x$table, 10), row.names = FALSE)
  invisible(x)
}

#' Infer adjacency matrix/matrices from a saved biSBM fit
#'
#' @param file        Path to the .rda file containing `bestSolutionAtQ`.
#' @param dataMatrix  The data matrix used in graph inference (n1 x n2 numeric matrix).
#' @param alpha       One or more alpha thresholds (numeric vector).
#' @param verbose     Logical; print summary info (default TRUE).
#'
#' @return If one alpha: an adjacency matrix.
#'         If multiple alphas: a named list of adjacency matrices.
#'
infer_adjacency_from_file <- function(file, dataMatrix, alpha = 0.10, best.id = NULL, verbose = TRUE) {
  stopifnot(file.exists(file))
  stopifnot(is.matrix(dataMatrix))
  
  # Load bestSolutionAtQ safely into its own environment
  env <- new.env(parent = emptyenv())
  loaded <- load(file, envir = env)
  
  if (!"bestSolutionAtQ" %in% loaded) {
    stop("The specified file must contain the object `bestSolutionAtQ`.")
  }
  
  bestSolutionAtQ <- get("bestSolutionAtQ", envir = env)
  if (!is.list(bestSolutionAtQ) || length(bestSolutionAtQ) == 0)
    stop("bestSolutionAtQ must be a non-empty list of model solutions.")
  
  # Identify the best solution by ICL
  ICLs <- sapply(bestSolutionAtQ, function(el) el$sbmParam$ICL)
  if (is.null(best.id)) best.id <- which.max(ICLs)
  currentSolution <- bestSolutionAtQ[[best.id]]
  
  # Extract parameters for graph inference
  clustering_row <- currentSolution$clustering_row
  clustering_col <- currentSolution$clustering_col
  theta          <- currentSolution$theta
  
  if (verbose) {
    message(sprintf("Selected solution #%d with ICL = %.3f (Q1=%s, Q2=%s, mean(pi)=%.3f)",
                    best.id, currentSolution$sbmParam$ICL,
                    currentSolution$sbmParam$Q1, currentSolution$sbmParam$Q2,
                    mean(theta$pi)))
  }
  
  # Run inference for one or multiple alpha thresholds
  alpha <- as.numeric(alpha)
  out_list <- lapply(alpha, function(a) {
    g <- graphInferenceNobiSBM(
      dataMatrix          = dataMatrix,
      nodeClusteringRow   = clustering_row,
      nodeClusteringCol   = clustering_col,
      theta               = theta,
      alpha               = a
    )
    g$A
  })
  
  # Name and return results
  names(out_list) <- paste0("alpha=", format(alpha, trim = TRUE))
  if (length(alpha) == 1L) {
    return(out_list[[1L]])
  } else {
    return(out_list)
  }
}

