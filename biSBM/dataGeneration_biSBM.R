#' simulation of a graph according the noisy bipartite stochastic block model
#' side one refers to rows while side two refers to columns.
#'
#' @param n1 number of nodes in side one
#' @param n2 number of nodes in side two
#' @param theta model parameters of the noisy bipartite stochastic block model
#' \describe{
#'   \item{pi}{connectivity parameters, Q1_Q2-vector}
#'   \item{alpha1}{latent block proportions for side one, Q1-vector}
#'   \item{alpha2}{latent block proportions for side two, Q2-vector}
#'   \item{nu0}{parameters of the null distribution}
#'   \item{nu}{parameters of the alternative distribution}
#' }
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#'
#' @return a list with:
#' \describe{
#'   \item{dataMatrix}{simulated matrix from the noisy bipartite stochastic block model}
#'   \item{theta}{model parameters of the noisy bipartite stochastic block model}
#'   \item{latentZ}{underlying latent node memberships}
#'   \item{latentAdj}{underlying latent binary graph}
#' }
#' @export
#'
#' @examples
#' n1 <- 10
#' n2 <- 15
#' Q1 <- 2
#' Q2 <- 2
#' theta <- list(alpha1 = c(0.5, 0.5), 
#'               alpha2 = c(0.4, 0.6),
#'               nu0 = c(0,.1),
#'               nu = list(mean = matrix(2, Q1, Q2), 
#'                         sd = matrix(0.1, Q1, Q2)),  
#'               pi=matrix(0.01,Q1,Q2) )
#' diag(theta$pi) <- 0.3
#' obs <- rnbisbm(n1, n2, theta, modelFamily='Gauss')
#' obs
rnbisbm <-  function(n1, n2, theta, modelFamily='Gauss'){
  N <- n1*n2
  Q1 <- length(theta$alpha1)
  Q2 <- length(theta$alpha2)
  
  # latent variables
  Z1 <- sample(1:Q1, n1, replace=TRUE, prob=theta$alpha1)
  Z2 <- sample(1:Q2, n2, replace=TRUE, prob=theta$alpha2)
  
  # adjacency matrix
  A <- matrix(0, n1, n2)
  for (i in 1:n1){
    A[i,] <- stats::rbinom(n2, 1, theta$pi[Z1[i],Z2])
  }
  
  # noisy observations under the null
  if (modelFamily=='Gauss'){
    X <- stats::rnorm(N, theta$nu0[1], theta$nu0[2])
  }
  if (modelFamily=='Gamma'){
    X <- stats::rgamma(N, shape = theta$nu0[1], rate = theta$nu0[2])
  }
  X <- matrix(X,n1,n2)
  
  for (i in 1:n1){
    nonzeroind <- which(A[i,]!=0)
    L <- length(nonzeroind)
    if (L>=1){
      if (modelFamily=='Gauss'){
        ind_i <- convertGroupPair(Z1[i], Z2[nonzeroind], Q1, Q2)
        X[i, nonzeroind] <- stats::rnorm(L, theta$nu$mean[ind_i], theta$nu$sd[ind_i])
      }
      if (modelFamily=='Gamma'){
        ind_i <- convertGroupPair(Z1[i], Z2[nonzeroind], Q1, Q2)
        X[i, nonzeroind] <- stats::rgamma(L, shape = theta$nu$mean[ind_i], rate = theta$nu$sd[ind_i])
      }
    }
  }
  
  return(list(dataMatrix=X, theta=theta, latentZ=list(Z1,Z2), latentAdj=A))
}

#' Function to generate (Fisher transformed) correlation matrix from log normal data
rcor_from_mvrnorm <- function(ns,n1,n2,Sigma,method='spearman'){
  logX <- MASS::mvrnorm(n=ns, mu=rep(0,n1+n2), Sigma=Sigma)
  X <- exp(logX[,1:n1])
  Y <- logX[,-(1:n1)]
  X.count <- round(X) # round the numbers to the nearest integer
  X.count[which(X.count==0)] <- 0.5
  X.rel <- sweep(X.count,1,rowSums(X.count), FUN='/')
  X.clr <- t(apply(X.rel, 1, function(a) log(a) - mean(log(a)))) # use this to correlate with Y
  
  rho <- cor(X.clr,Y,method = method)
  
  dataMatrix <- sqrt(ns-3)*(log(1+rho) - log(1-rho))/2
  return(dataMatrix)
  
}

#' Rarefy count matrix to equal total depths across samples
#'
#' @param X A numeric matrix of raw counts (rows = samples, columns = taxa/features).
#' @return A matrix with counts rescaled so that all samples have the same total depth
#'         equal to the minimum observed sequencing depth (rounded to integers).
#' @examples
#' set.seed(1)
#' X <- matrix(rpois(20, lambda = 10), nrow = 4)
#' X.rarefied <- rarefy_to_min_depth(X)
#' rowSums(X.rarefied)
rarefy_to_min_depth <- function(X) {
  if (!is.matrix(X))
    stop("Input X must be a numeric matrix with samples in rows and features in columns.")
  if (any(X < 0))
    stop("Counts must be nonnegative.")
  
  depths <- rowSums(X)
  min_depth <- min(depths)
  
  if (min_depth == 0)
    stop("At least one sample has zero total counts; cannot normalize.")
  
  # Normalize counts to relative abundance, rescale to min depth, and round
  X.norm <- sweep(X, 1, depths, FUN = "/")
  X.equalized <- round(X.norm * min_depth)
  
  return(X.equalized)
}


#' Mixed-marginal synthetic data from joint Gaussian copula
#'
#' Generate synthetic microbiome counts (NB) and metabolite abundances (log-normal)
#' using a single joint correlation matrix (Sigma). NB marginals are estimated
#' from X (taxon counts); log-normal marginals are estimated from Y (metabolites).
#'
#' @param X  (n x p_x) matrix/data.frame of taxon counts used to fit NB marginals.
#' @param Y  (n x p_y) matrix/data.frame of metabolite abundances used to fit LN marginals.
#' @param Sigma  ((p_x + p_y) x (p_x + p_y)) target latent Gaussian correlation matrix.
#'               The first p_x dimensions correspond to taxa (X), the remaining p_y to metabolites (Y).
#' @param n  number of samples to simulate (default: nrow(X) if available, else nrow(Y)).
#' @param eps small positive number added when taking logs for LN fitting if zeros present (default 1e-6).
#' @param cap_large_size large "size" for NB when var <= mean (Poisson-like; default 1e6).
#' @param return_latent if TRUE, also return the latent Gaussian draws Z and uniforms U (default FALSE).
#' @param seed optional RNG seed for reproducibility.
#'
#' @return list with elements:
#'   - X: (n x p_x) synthetic taxon counts (NB marginals)
#'   - Y: (n x p_y) synthetic metabolite abundances (log-normal marginals)
#'   - params_nb: data.frame of NB (mu, size) used per taxon
#'   - params_ln: data.frame of LN (meanlog, sdlog) used per metabolite
#'   - (optional) Z, U if return_latent = TRUE
#'
#' @details
#' Uses a Gaussian copula / NorTA approach:
#'   1) Draw Z ~ N(0, Sigma).
#'   2) U = Φ(Z) to get uniform[0,1] margins.
#'   3) For taxa: X_ij = qnbinom(U_ij; size_j, mu_j) with NB params fit from column j of X.
#'   4) For metabolites: Y_ik = qlnorm(U_i, meanlog_k, sdlog_k) with LN params fit from column k of Y.
#'
#' NB fitting is method-of-moments: size = mu^2 / (var - mu), capped at cap_large_size if var <= mu.
#' LN fitting uses sample mean/SD on log scale; if non-positive values exist, adds eps before log.
#'
#' @importFrom MASS mvrnorm
#' @examples
#' # Suppose X (counts), Y (metabolites) already exist, and Sigma is (ncol(X)+ncol(Y))^2:
#' # sim <- synth_comm_from_XY(X, Y, Sigma)
synth_comm_from_XY <- function(X, Y, Sigma, n,
                               eps = 1e-6,
                               library_size = NULL,
                               library_scale = rep(1,n),
                               cap_large_size = 1e6,
                               return_latent = FALSE,
                               seed = NULL) {
  # -- Coerce to matrices
  X <- as.matrix(X)
  depths <- rowSums(X) # raw counts data, rows are obs, n by p
  X.n  <- sweep(X, 1, depths, FUN = '/')
  if (is.null(library_size)) library_size = min(depths)
  X <- round(X.n * library_size)
  Y <- as.matrix(Y)
  
  p_x <- ncol(X)
  p_y <- ncol(Y)
  if (is.null(p_x)) p_x <- 0
  if (is.null(p_y)) p_y <- 0
  p <- p_x + p_y
  
  if (!is.matrix(Sigma) || any(dim(Sigma) != c(p, p))) {
    stop("Sigma must be a (p_x + p_y) x (p_x + p_y) matrix matching cbind(X, Y) columns.")
  }
  
  # -- Basic checks
  if (p == 0) stop("No variables: ncol(X) + ncol(Y) must be > 0.")
  if (!is.numeric(Sigma) || !all(is.finite(Sigma))) stop("Sigma must be numeric and finite.")
  # Symmetrize tiny asymmetries and try to coerce to correlation
  # (Assume user passes a proper correlation; we lightly guard.)
  Sigma <- (Sigma + t(Sigma)) / 2
  d <- sqrt(diag(Sigma))
  if (any(d <= 0)) stop("Sigma must have positive diagonal.")
  Sigma <- Sigma / (d %o% d) # scale to correlation
  
  # -- Fit NB marginals from X
  fit_nb_one <- function(x, cap_large_size) {
    mu <- mean(x, na.rm = TRUE)
    v  <- stats::var(x, na.rm = TRUE)
    if (!is.finite(mu) || length(x) <= 1) {
      return(c(mu = 0, size = cap_large_size))
    }
    if (!is.finite(v)) v <- 0
    if (mu <= 0) {
      # All zeros (or nonpositive); keep degenerate at zero via huge size and mu=0
      return(c(mu = 0, size = cap_large_size))
    }
    if (v > mu) {
      size <- mu^2 / (v - mu)
      if (!is.finite(size) || size <= 0) size <- cap_large_size
    } else {
      # Poisson or underdispersed: approximate with large size
      size <- cap_large_size
    }
    c(mu = mu, size = size)
  }
  
  params_nb <- NULL
  if (p_x > 0) {
    nb_mat <- t(apply(X, 2, fit_nb_one, cap_large_size = cap_large_size))
    params_nb <- data.frame(feature = colnames(X) %||% paste0("taxon_", seq_len(p_x)),
                            mu = nb_mat[, "mu"], size = nb_mat[, "size"],
                            row.names = NULL, check.names = FALSE)
  }
  
  # -- Fit LN marginals from Y
  fit_ln_one <- function(y, eps) {
    y_adj <- as.numeric(y)
    # If any non-positive, add eps before log to avoid -Inf
    if (any(y_adj <= 0 | !is.finite(y_adj))) {
      y_adj[!is.finite(y_adj)] <- NA_real_
      y_adj <- y_adj + eps
    }
    ly <- log(y_adj)
    ly <- ly[is.finite(ly)]
    if (length(ly) < 2) {
      # fallback: meanlog ~ log(mean+eps), sdlog small
      m <- log(mean(y_adj, na.rm = TRUE) + eps)
      s <- 1e-6
    } else {
      m <- mean(ly, na.rm = TRUE)
      s <- stats::sd(ly, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) s <- 1e-6
    }
    c(meanlog = m, sdlog = s)
  }
  
  params_ln <- NULL
  if (p_y > 0) {
    ln_mat <- t(apply(Y, 2, fit_ln_one, eps = eps))
    params_ln <- data.frame(feature = colnames(Y) %||% paste0("metab_", seq_len(p_y)),
                            meanlog = ln_mat[, "meanlog"], sdlog = ln_mat[, "sdlog"],
                            row.names = NULL, check.names = FALSE)
  }
  
  # -- Draw latent Gaussian and transform
  if (!is.null(seed)) set.seed(seed)
  # numeric stability: small jitter to ensure PD if needed
  # (user should pass PD Sigma; we avoid nearPD dependency here)
  Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma, empirical = FALSE)
  U <- pnorm(Z)
  
  # -- Apply marginals
  X_sim <- NULL
  if (p_x > 0) {
    X_sim <- matrix(0L, nrow = n, ncol = p_x)
    for (j in seq_len(p_x)) {
      mu_j   <- params_nb$mu[j]
      size_j <- params_nb$size[j]
      if (mu_j <= 0) {
        X_sim[, j] <- 0L
      } else {
        # qnbinom uses 'size' and 'mu'
        X_sim[, j] <- stats::qnbinom(U[, j], size = size_j, mu = mu_j)
      }
    }
    colnames(X_sim) <- params_nb$feature
    storage.mode(X_sim) <- "integer"
  }
  X_sim <- sweep(X_sim, 1, library_scale, FUN="*")
  
  Y_sim <- NULL
  if (p_y > 0) {
    Y_sim <- matrix(NA_real_, nrow = n, ncol = p_y)
    for (k in seq_len(p_y)) {
      meanlog_k <- params_ln$meanlog[k]
      sdlog_k   <- params_ln$sdlog[k]
      Y_sim[, k] <- stats::qlnorm(U[, p_x + k], meanlog = meanlog_k, sdlog = sdlog_k)
    }
    colnames(Y_sim) <- params_ln$feature
  }
  
  out <- list(X = X_sim, Y = Y_sim, params_nb = params_nb, params_ln = params_ln)
  if (return_latent) out[c("Z", "U")] <- list(Z = Z, U = U)
  out
}

# helper: null-coalescing for names
`%||%` <- function(a, b) if (!is.null(a)) a else b
t_to_z <- function(tstat, df) {
  # two-sided
  p <- 2 * pt(-abs(tstat), df)
  z <- sign(tstat) * qnorm(1 - p / 2)
  return(z)
}

#' Function to generate multivariate microbiome--metabolomic data from a lognormal distribution. 
#' In addition, we add a count sampling layer for the microbiome data. 
#' @param n sample size
#' @param p1 number of microbes
#' @param mu mean of the latent basis
#' @param Sigma covariance of the latent basis, the first p1 are microbes.
#' @param library_scale a vector of length n that provides the library size for each sample. 
#' @return A list with the latent basis (`basis`), observed counts (`taxa`), and observed relative abundances (`rel`)
para_data_generation <- function(n, p1, mu, Sigma, library_scale){
  Sigma_mat_decomp = t(chol(Sigma))
  p <- nrow(Sigma) # total dimensions
  x = Sigma_mat_decomp %*% t(matrix(rnorm(n=n*p, mean = 0), n, p, byrow=T)) + mu
  y = t(exp(x)) # log normal data
  
  # add a count sampling layer to the microbiome data
  comp_mat = sweep(y[,1:p1],1,STATS = rowSums(y[,1:p1]), FUN='/')
  data = t(
    do.call(
      cbind,
      lapply(1:n, 
             function(i) rmultinom(n=1,size=library_scale[i], prob = comp_mat[i,]))
    )
  )
  obs_rel = sweep(data,1,STATS = rowSums(y), FUN='/')
  
  return(list(Z = t(x),
              rel = obs_rel, 
              X = data,
              Y = y[,-(1:p1)]))
}

para_data_generation_ln <- function(n, p1, mu, Sigma){
  Sigma_mat_decomp = t(chol(Sigma))
  p <- nrow(Sigma) # total dimensions
  x = Sigma_mat_decomp %*% t(matrix(rnorm(n=n*p, mean = 0), n, p, byrow=T)) + mu
  y = t(exp(x)) # log normal data
  counts <- round(y[,1:p1])
  
  # add a count sampling layer to the microbiome data
  comp_mat = sweep(counts,1,STATS = rowSums(counts), FUN='/')
  data = comp_mat
  obs_rel = sweep(data,1,STATS = rowSums(data), FUN='/')
  
  return(list(Z = t(x),
              rel = obs_rel, 
              X = data,
              counts = counts,
              Y = y[,-(1:p1)]))
}

# --------------------------------------------------------------------
# Joint correlation from bipartite network with data-driven diagonals
# --------------------------------------------------------------------

#' Build a joint correlation matrix for (X, Y)
#' - R_xx inferred from data block X
#' - R_yy inferred from data block Y
#' - R_xy from bipartite network A, optionally modulated by cluster pairs (Zx, Zy)
#'
#' @param A      px x py bipartite adjacency (0/1 or weights)
#' @param Zx     integer cluster labels for X variables (length px)
#' @param Zy     integer cluster labels for Y variables (length py)
#' @param X      n x px data matrix (used to infer R_xx). If NULL, falls back to identity/random.
#' @param Y      n x py data matrix (used to infer R_yy). If NULL, falls back to identity/random.
#' @param cross_rho  base correlation magnitude for edges in A (default 0.3)
#' @param rho_by_block  optional kx x ky matrix with per-(xblock,yblock) multipliers
#'                      applied to cross_rho (i.e., edge rho = cross_rho * rho_by_block[Zx,Zy]).
#'                      If NULL, uses 1 for every active block pair.
#' @param method_within  how to infer within-block correlations: "data", "shrink", "identity", "random"
#' @param corr_type      "pearson" or "spearman" when using data
#' @param shrink_lambda  if method_within="shrink", ridge shrinkage level to I (0..1)
#' @param random_eta     strength for random corr (smaller => stronger off-diagonals)
#' @param jitter_xy      small N(0, sd=jitter_xy) noise inside active cross edges
#' @param project_nearPD logical; if TRUE, final nearPD(corr=TRUE)
#'
#' @return list with R (joint correlation), Rxx, Ryy, Rxy (pre-projection), scale_factor, eigvals_before
#'
build_corr_from_bipartite <- function(
    A, Zx, Zy,
    X = NULL, Y = NULL,
    cross_rho = 0.3,
    rho_by_block = NULL,
    method_within = c("data", "shrink", "identity", "random"),
    corr_type = c("pearson", "spearman"),
    shrink_lambda = 0.1,
    random_eta = 1.5,
    jitter_xy = 0.00,
    project_nearPD = TRUE,
    # NEW:
    cross_mode = c("scale", "whitened", "factor"),
    target_spec = 0.98,      # for "whitened" / "factor"
    factor_rank = 1          # for "factor": number of shared factors
) {
  method_within <- match.arg(method_within)
  corr_type <- match.arg(corr_type)
  cross_mode <- match.arg(cross_mode)
  
  A <- as.matrix(A)
  px <- nrow(A); py <- ncol(A)
  stopifnot(length(Zx) == px, length(Zy) == py)
  
  # --- helpers (same as before) ---
  .random_corr <- function(p, eta = 1.5) {
    if (requireNamespace("ClusterGeneration", quietly = TRUE)) {
      ClusterGeneration::rcorrmatrix(p, alphad = eta)
    } else { M <- matrix(rnorm(p*p), p,p); S <- crossprod(M); D <- sqrt(diag(S)); S/(D%o%D) }
  }
  make_within <- function(dat, p) {
    if (method_within %in% c("data","shrink")) {
      if (is.null(dat)) stop("X/Y required when method_within='", method_within, "'.")
      Rhat <- suppressWarnings(stats::cor(dat, method = corr_type, use = "pairwise.complete.obs"))
      Rhat[!is.finite(Rhat)] <- 0; diag(Rhat) <- 1
      if (method_within == "shrink") (1 - shrink_lambda) * Rhat + shrink_lambda * diag(p) else Rhat
    } else if (method_within == "identity") diag(p) else .random_corr(p, eta = random_eta)
  }
  # Build a raw cross-block from A (+clusters)
  make_raw_cross <- function() {
    Rxy <- matrix(0, px, py)
    if (is.null(rho_by_block)) rho_by_block <- matrix(1, nrow = max(Zx), ncol = max(Zy))
    for (i in seq_len(px)) for (j in seq_len(py)) if (!is.na(A[i,j]) && A[i,j] != 0) {
      val <- cross_rho * rho_by_block[Zx[i], Zy[j]]
      if (jitter_xy > 0) val <- val + rnorm(1, sd = jitter_xy)
      val <- val * as.numeric(A[i,j])
      Rxy[i,j] <- max(min(val, 0.99), -0.99)
    }
    Rxy
  }
  eig_pd <- function(M) {
    ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA)
    if (any(!is.finite(ev)) || min(ev) <= 1e-8) M + 1e-6*diag(nrow(M)) else M
  }
  inv_sqrt <- function(M) { ee <- eigen(M, symmetric=TRUE); v <- pmax(ee$values,1e-8)
  ee$vectors %*% diag(1/sqrt(v), length(v)) %*% t(ee$vectors) }
  sqrtm <- function(M) { ee <- eigen(M, symmetric=TRUE); v <- pmax(ee$values,0)
  ee$vectors %*% diag(sqrt(v), length(v)) %*% t(ee$vectors) }
  
  # 1) Within-blocks
  Rxx <- eig_pd(make_within(X, px))
  Ryy <- eig_pd(make_within(Y, py))
  
  # 2) Cross block
  if (cross_mode == "scale") {
    Rxy_raw <- make_raw_cross()
    Rxm <- inv_sqrt(Rxx); Rym <- inv_sqrt(Ryy)
    Tmat <- Rxm %*% Rxy_raw %*% Rym
    smax <- tryCatch(base::svd(Tmat, nu=0, nv=0)$d[1], error=function(e) 0)
    scale_factor <- if (smax > 0) min(1, 0.97/smax) else 1
    Rxy <- scale_factor * Rxy_raw
  } else if (cross_mode == "whitened") {
    # Build T in whitened space, normalize to target_spec, then unwhiten
    Rxy_raw <- make_raw_cross()
    Rxm <- inv_sqrt(Rxx); Rym <- inv_sqrt(Ryy)
    T_base <- Rxm %*% Rxy_raw %*% Rym
    smax <- tryCatch(base::svd(T_base, nu=0, nv=0)$d[1], error=function(e) 0)
    if (smax == 0) {
      T <- T_base
    } else {
      T <- (target_spec / smax) * T_base  # directly hit desired spectral norm
    }
    Sx <- sqrtm(Rxx); Sy <- sqrtm(Ryy)
    Rxy <- Sx %*% T %*% Sy
    scale_factor <- target_spec / (smax + (smax==0))
  } else { # cross_mode == "factor"
    # Low-rank coupling: Rxy = Sx U diag(g) V' Sy, with singular values g < 1
    Sx <- sqrtm(Rxx); Sy <- sqrtm(Ryy)
    # Build block-constrained loadings U,V (normalized to blocks implied by A)
    gx_ids <- Zx; gy_ids <- Zy
    kx <- max(gx_ids); ky <- max(gy_ids)
    # Make U,V by averaging within active blocks from A
    U <- matrix(0, px, factor_rank); V <- matrix(0, py, factor_rank)
    # simple scheme: first factor loads on the densest active pairs
    degx <- rowSums(A != 0); degy <- colSums(A != 0)
    U[,1] <- (degx + 1e-8); U[,1] <- U[,1]/sqrt(sum(U[,1]^2))
    V[,1] <- (degy + 1e-8); V[,1] <- V[,1]/sqrt(sum(V[,1]^2))
    if (factor_rank > 1) {  # crude orthogonal fillers
      set.seed(1)
      U[,2:factor_rank] <- qr.Q(qr(cbind(U[,1], matrix(rnorm(px*(factor_rank-1)), px))), complete=FALSE)[,2:factor_rank, drop=FALSE]
      V[,2:factor_rank] <- qr.Q(qr(cbind(V[,1], matrix(rnorm(py*(factor_rank-1)), py))), complete=FALSE)[,2:factor_rank, drop=FALSE]
    }
    g <- rep(target_spec / factor_rank, factor_rank)  # keep spectral norm < 1
    Rxy <- Sx %*% (U %*% diag(g, factor_rank) %*% t(V)) %*% Sy
    scale_factor <- target_spec
    # Optionally zero out entries where A==0 to respect sparsity, then nearPD later:
    Rxy[A == 0] <- 0
  }
  
  # 3) Assemble + optional nearPD
  R <- matrix(0, px+py, px+py)
  R[1:px,1:px] <- Rxx
  R[(px+1):(px+py),(px+1):(px+py)] <- Ryy
  R[1:px,(px+1):(px+py)] <- Rxy
  R[(px+1):(px+py),1:px] <- t(Rxy)
  R <- (R + t(R))/2; diag(R) <- 1
  
  eigvals_before <- tryCatch(eigen(R, symmetric=TRUE, only.values=TRUE)$values, error=function(e) NA)
  if (project_nearPD && requireNamespace("Matrix", quietly = TRUE)) {
    R <- as.matrix(Matrix::nearPD(R, corr=TRUE)$mat)
  }
  
  list(R=R, Rxx=Rxx, Ryy=Ryy, Rxy=Rxy, mode=cross_mode,
       scale_factor=scale_factor, eigvals_before=eigvals_before)
}


generate_covariance <- function(A, Z, n1, n2, safety = 1-1e-03, signed = TRUE){
  # Ensure PSD of the full covariance by shrinking the spectral norm of Sxy
  var_x <- 1
  var_y <- 1
  # safety <- 
  s2_max <- sqrt(var_x * var_y)
  if (signed){
    Q1 <- length(unique(Z[[1]]))
    Q2 <- length(unique(Z[[2]]))
    for (j in 1:Q1){
      for (k in 1:Q2){
        A[Z[[1]]==j, Z[[2]]==k] <- A[Z[[1]]==j, Z[[2]]==k] * (2*rbinom(1,1,0.5) - 1)
      }
    }
  }
  Sxy <- A
  sp <- if (all(Sxy == 0)) 0 else svd(Sxy, nu = 0, nv = 0)$d[1]
  if (sp >= s2_max * safety) {
    Sxy <- Sxy * ((s2_max * safety) / (sp + 1e-12))
  }
  
  # Assemble covariance with diagonal within-blocks
  Sxx <- diag(var_x, n1)
  Syy <- diag(var_y, n2)
  rbind(cbind(Sxx, Sxy),
        cbind(t(Sxy), Syy))
  
}
