# stability selection

# ---- helpers ---------------------------------------------------------------
.argmax <- function(x) {
  # x is a numeric vector or a matrix
  if (is.matrix(x)) max.col(t(x)) else which.max(x)
}

# robust mclapply: falls back to lapply on Windows
.mapply_parallel <- function(X, FUN, mc.cores = 1, ...) {
  if (.Platform$OS.type == "windows" || mc.cores <= 1) {
    lapply(X, FUN, ...)
  } else {
    parallel::mclapply(X, FUN, mc.cores = mc.cores, ...)
  }
}

# extract hard clusters from a fit list element
.get_hard_labels <- function(fit_el) {
  cr <- fit_el$clustering_row
  cc <- fit_el$clustering_col
  # If not present, derive from variational params:
  if (is.null(cr) && !is.null(fit_el$sbmParam$clusterProba[[1]])) {
    cr <- apply(fit_el$sbmParam$clusterProba[[1]], 2, .argmax)
  }
  if (is.null(cc) && !is.null(fit_el$sbmParam$clusterProba[[2]])) {
    cc <- apply(fit_el$sbmParam$clusterProba[[2]], 2, .argmax)
  }
  list(row = cr, col = cc,
       Q1 = fit_el$sbmParam$Q1, Q2 = fit_el$sbmParam$Q2,
       ICL = fit_el$sbmParam$ICL)
}

# choose a single "best" solution among the list returned by fitNobiSBM
.pick_solution <- function(fit_list, select_by = c("ICL","first"), Q1 = NULL, Q2 = NULL) {
  select_by <- match.arg(select_by)
  if (!is.null(Q1) && !is.null(Q2)) {
    idx <- which(sapply(fit_list, function(el) el$sbmParam$Q1 == Q1 &&
                          el$sbmParam$Q2 == Q2))
    if (length(idx)) return(fit_list[[idx[1]]])
  }
  if (select_by == "ICL") {
    icl <- sapply(fit_list, function(el) el$sbmParam$ICL)
    return(fit_list[[which.max(icl)]])
  } else {
    return(fit_list[[1]])
  }
}

#' Update Consensus Counts (Internal)
#'
#' @description
#' Internal helper to update numerator/denominator consensus matrices with one set
#' of indices and labels. Returns the updated matrices (R is pass-by-value).
#'
#' @param cons_same Integer matrix; numerator counts for same-cluster pairs.
#' @param cons_seen Integer matrix; denominator counts for co-observed pairs.
#' @param idx Integer vector of in-subsample indices (relative to full data).
#' @param labels Integer labels for \code{idx}.
#'
#' @return List with updated \code{cons_same} and \code{cons_seen}.
#' @keywords internal
#' @noRd
.update_consensus <- function(cons_same, cons_seen, idx, labels) {
  if (length(idx) < 2L) return(list(cons_same = cons_same, cons_seen = cons_seen))
  
  # co-observed pairs: +1 to all pairs in the submatrix
  cons_seen[idx, idx] <- cons_seen[idx, idx] + 1L
  
  # same-cluster pairs: +1 within each cluster block
  ord <- order(labels)
  idx  <- idx[ord]
  lab  <- labels[ord]
  rle_lab <- rle(lab)
  start <- 1L
  for (k in seq_along(rle_lab$lengths)) {
    L <- rle_lab$lengths[k]
    if (L >= 2L) {
      block <- idx[start:(start + L - 1L)]
      cons_same[block, block] <- cons_same[block, block] + 1L
    }
    start <- start + L
  }
  list(cons_same = cons_same, cons_seen = cons_seen)
}

# summarize stability for a given partition using a consensus matrix
.partition_stability <- function(consensus, labels) {
  # consensus: in [0,1]; labels length = n
  n <- length(labels)
  if (n <= 1) return(list(overall_mean = NA_real_, overall_median = NA_real_,
                          by_cluster = data.frame(cluster=integer(), mean=double(), median=double())))
  # within-cluster consensus (exclude diagonal)
  res <- lapply(split(seq_len(n), labels), function(idx) {
    if (length(idx) <= 1) return(c(mean = NA_real_, median = NA_real_))
    M <- consensus[idx, idx, drop = FALSE]
    v <- M[upper.tri(M, diag = FALSE)]
    c(mean = mean(v), median = median(v))
  })
  byc <- do.call(rbind, res)
  data.frame(cluster = as.integer(rownames(byc)), mean = byc[,1], median = byc[,2], row.names = NULL)
}

# ---- metrics on one consensus + labels -------------------------------------
#' Proportion of Ambiguous Clustering (PAC)
#'
#' @description
#' Computes the fraction of off-diagonal consensus entries in an ambiguity band
#' \code{(lower, upper)}. Lower PAC indicates crisper clustering.
#'
#' @param M Square consensus matrix in \eqn{[0,1]}.
#' @param lower,upper Numeric scalars defining the ambiguity band; defaults \code{0.1}, \code{0.9}.
#'
#' @return Numeric scalar in \eqn{[0,1]}.
#'
#' @examples
#' \dontrun{
#' pac <- consensus_PAC(stab$consensus_row)
#' }
#' @export
consensus_PAC <- function(M, lower = 0.1, upper = 0.9) {
  # proportion of off-diagonal entries in (lower, upper)
  stopifnot(nrow(M) == ncol(M))
  v <- M[upper.tri(M, diag = FALSE)]
  mean(v > lower & v < upper, na.rm = TRUE)
}

#' Within- and Between-Cluster Consensus Summaries
#'
#' @description
#' Returns mean and median of consensus values within clusters and between clusters,
#' given hard labels.
#'
#' @param M Square consensus matrix in \eqn{[0,1]}.
#' @param labels Integer or factor vector of cluster assignments.
#'
#' @return Named numeric vector with elements:
#' \code{within_mean}, \code{within_median}, \code{between_mean}, \code{between_median}.
#'
#' @note Requires at least two clusters to have a defined between-cluster summary.
#'
#' @examples
#' \dontrun{
#' wb_row <- consensus_within_between(stab$consensus_row, stab$full_fit$clustering_row)
#' }
#' @export
consensus_within_between <- function(M, labels) {
  n <- length(labels)
  stopifnot(nrow(M) == n, ncol(M) == n)
  # within-cluster
  W <- unlist(lapply(split(seq_len(n), labels), function(idx) {
    if (length(idx) <= 1) return(NA_real_)
    m <- M[idx, idx, drop = FALSE]
    m[upper.tri(m, diag = FALSE)]
  }))
  # between-cluster
  B <- unlist({
    labs <- unique(labels)
    out <- list()
    if (length(labs) >= 2) {
      for (i in 1:(length(labs)-1)) for (j in (i+1):length(labs)) {
        I <- which(labels == labs[i]); J <- which(labels == labs[j])
        out[[length(out)+1]] <- as.vector(M[I, J, drop = FALSE])
      }
    }
    out
  })
  c(within_mean = mean(W, na.rm = TRUE),
    within_median = stats::median(W, na.rm = TRUE),
    between_mean = mean(B, na.rm = TRUE),
    between_median = stats::median(B, na.rm = TRUE))
}

#' Silhouette on Consensus Distances
#'
#' @description
#' Computes the average silhouette width using \code{1 - consensus} as a distance.
#'
#' @param M Square consensus matrix in \eqn{[0,1]}.
#' @param labels Integer or factor vector of cluster assignments.
#' @param method Linkage method used inside \code{hclust} if needed; kept for symmetry.
#'
#' @return Average silhouette width (numeric). Returns \code{NA} if only one cluster
#' or if \pkg{cluster} is not installed.
#'
#' @examples
#' \dontrun{
#' sil_row <- consensus_silhouette(stab$consensus_row, stab$full_fit$clustering_row)
#' }
#' @importFrom stats as.dist
#' @export
consensus_silhouette <- function(M, labels, method = "average") {
  if (!requireNamespace("cluster", quietly = TRUE))
    return(NA_real_)
  # distance = 1 - consensus; clamp for safety
  D <- 1 - pmin(pmax(M, 0), 1)
  d <- stats::as.dist(D)
  sil <- cluster::silhouette(as.integer(labels), d)
  mean(sil[, "sil_width"])
}


#' Composite Stability Score (Consensus + Geometry; bigger is better)
#'
#' @description
#' Combines multiple stability indicators—within/between consensus, silhouette, and PAC—
#' into a single score. Higher is better. Handles the single-cluster case gracefully.
#'
#' @param M Square consensus matrix in \eqn{[0,1]}.
#' @param labels Hard labels. If only one cluster is present, silhouette is undefined;
#' the score reduces to \code{within * (1 - PAC)}.
#' @param w_within,w_between,w_sil,w_pac Weights for components (sum need not be 1).
#' @param pac_lower,pac_upper Band for PAC ambiguity.
#'
#' @return List with elements \code{score}, \code{within_mean}, \code{between_mean},
#' \code{silhouette}, \code{PAC}.
#'
#' @examples
#' \dontrun{
#' sc_row <- stability_score(stab$consensus_row, stab$full_fit$clustering_row)
#' sc_col <- stability_score(stab$consensus_col, stab$full_fit$clustering_col)
#' mean(c(sc_row$score, sc_col$score))
#' }
#' @seealso \code{\link{consensus_PAC}}, \code{\link{consensus_silhouette}}
#' @export
stability_score <- function(M, labels,
                            w_within = 0.4, w_between = 0.2,
                            w_sil = 0.3, w_pac = 0.1,
                            pac_lower = 0.1, pac_upper = 0.9) {
  # if only one cluster, skip metrics that require >= 2 clusters
  n_clusters <- length(unique(labels))
  if (n_clusters < 2) {
    within <- mean(M[upper.tri(M)], na.rm = TRUE)
    pac <- consensus_PAC(M, pac_lower, pac_upper)
    # Assign a neutral score: stability = within * (1 - pac)
    score <- within * (1 - pac)
    return(list(score = score,
                within_mean = within,
                between_mean = NA,
                silhouette = NA,
                PAC = pac))
  }
  
  # ---- regular case (>= 2 clusters) ----
  wb <- consensus_within_between(M, labels)
  sil <- consensus_silhouette(M, labels)
  pac <- consensus_PAC(M, pac_lower, pac_upper)
  
  within <- wb["within_mean"]
  between <- 1 - wb["between_mean"]
  sil_ok <- ifelse(is.na(sil), 0, pmin(pmax((sil + 1)/2, 0), 1))
  pac_ok <- 1 - pac
  
  score <- w_within*within + w_between*between + w_sil*sil_ok + w_pac*pac_ok
  
  list(score = as.numeric(score),
       within_mean = within,
       between_mean = 1 - between,
       silhouette = sil,
       PAC = pac)
}

# ---- main wrapper -----------------------------------------------------------
#' Stability via Subsampling for Noisy Bipartite SBM
#'
#' @description
#' Runs stability analysis for a noisy bipartite stochastic block model (nobiSBM)
#' by repeatedly subsampling rows and columns, refitting the model, and building
#' row/column **consensus co-clustering** matrices. Optionally shows a progress bar
#' (serial or parallel) using \pkg{pbapply}.
#'
#' @param dataMatrix Numeric matrix, \eqn{n_1 \times n_2}. Required.
#' @param model Character scalar. Model family passed to \code{fitNobiSBM()} (e.g., "Gauss0", "Gauss01", "Exp1").
#' @param sbmSize List with elements \code{Q1}, \code{Q2} (vectors of candidate block counts)
#' and \code{explor} (numeric, exploration multiplier). If \code{fix_Q1}/\code{fix_Q2} are set,
#' they take precedence over the grid.
#' @param exclusive_Rows_Cols Logical; if \code{TRUE}, the numbers of row and column blocks are set to be equal.
#' @param num.iter Integer; maximum iterations per fit.
#' @param initParam List of initialization controls passed to \code{fitNobiSBM()}.
#' @param B Integer; number of subsamples (bootstrap iterations). Default \code{50}.
#' @param row_frac,col_frac Numeric in \code{(0,1]}; fraction of rows/cols drawn without replacement per subsample.
#' @param select_by \code{"ICL"} or \code{"first"}; how to pick the best solution from the list returned by \code{fitNobiSBM()}.
#' @param fix_Q1,fix_Q2 Optional integers. If set, fix the number of row/column blocks for every subsample.
#' @param seed Optional integer for reproducibility.
#' @param keep_per_run Logical; if \code{TRUE}, stores per-subsample indices/labels and chosen \code{(Q1,Q2,ICL)}.
#' @param verbose Logical; if \code{TRUE}, shows a progress bar (requires \pkg{pbapply} for parallel progress).
#'
#' @details
#' The function builds two consensus matrices:
#' \itemize{
#'   \item \code{consensus_row[i,j]} = proportion of subsamples where rows \code{i} and \code{j}
#'         are both present and assigned to the same row-block.
#'   \item \code{consensus_col[i,j]} = analogous for columns.
#' }
#' Diagonals are set to 1. Pairs never co-observed keep a 0 denominator and are left at 0
#' in the final consensus (recorded in \code{seen_*} for auditing).
#'
#' A single full-data fit is run at the end (respecting \code{fix_Q1}/\code{fix_Q2} if given)
#' and used to summarize per-cluster stability.
#'
#' @return
#' An object of class \code{"nobiSBM_stability"} with elements:
#' \itemize{
#'   \item \code{consensus_row}, \code{consensus_col} (\eqn{[0,1]} symmetric matrices)
#'   \item \code{seen_row}, \code{seen_col} (integer matrices of co-observation counts)
#'   \item \code{stability_summary}: overall and per-cluster mean/median consensus
#'   \item \code{full_fit}: best full-data fit (as returned by \code{fitNobiSBM()})
#'   \item \code{settings}: list of key arguments used
#'   \item \code{per_run}: list of per-iteration results (present if \code{keep_per_run=TRUE})
#' }
#'
#' @section Edge cases:
#' If \code{Q1=1} or \code{Q2=1}, some stability metrics (e.g., silhouette) are undefined;
#' see \code{\link{stability_score}} for recommended handling.
#'
#' @examples
#' \dontrun{
#' stab <- stability_nobiSBM(
#'   dataMatrix = d$dataMatrix,
#'   model = "Gauss01",
#'   sbmSize = list(Q1 = 2:5, Q2 = 2:5, explor = 1.5),
#'   B = 50, row_frac = 0.8, col_frac = 0.8,
#'   nbCores = 4, seed = 1, verbose = TRUE
#' )
#' }
#'
#' @seealso \code{\link{plot_consensus_heatmap}}, \code{\link{stability_score}},
#' \code{\link{score_stability_grid}}
#' @importFrom stats as.dist hclust cutree median
#' @export
stability_nobiSBM <- function(
    dataMatrix,
    model = 'Gauss0',
    sbmSize = list(Q1 = 1:5, Q2 = 1:5, explor = 1.5),
    exclusive_Rows_Cols = FALSE,
    num.iter = 5,
    initParam = list(nbOfbeta = NULL, nbOfPointsPerbeta = NULL,
                     maxNbOfPasses = NULL, minNbOfPasses = 1),
    B = 50,
    row_frac = 0.8,
    col_frac = 0.8,
    select_by = c("ICL","first"),
    fix_Q1 = NULL,
    fix_Q2 = NULL,
    seed = NULL,
    keep_per_run = FALSE,
    verbose = TRUE
) {
  select_by <- match.arg(select_by)
  if (!is.null(seed)) set.seed(seed)
  
  n1 <- nrow(dataMatrix); n2 <- ncol(dataMatrix)
  if (row_frac <= 0 || row_frac > 1 || col_frac <= 0 || col_frac > 1)
    stop("row_frac and col_frac must be in (0,1].")
  m1 <- max(1L, floor(row_frac * n1))
  m2 <- max(1L, floor(col_frac * n2))
  
  cons_same_row <- matrix(0L, n1, n1)
  cons_seen_row <- matrix(0L, n1, n1)
  cons_same_col <- matrix(0L, n2, n2)
  cons_seen_col <- matrix(0L, n2, n2)
  per_run <- if (keep_per_run) vector("list", B) else NULL
  
  runs <- seq_len(B)
  
  # ---------- progress + parallel via pbapply ----------
  use_pb <- verbose && requireNamespace("pbapply", quietly = TRUE)
  
  # Worker function: returns only what we need to aggregate later
  worker <- function(b, dataMatrix, n1, n2, m1, m2,
                     model, sbmSize, exclusive_Rows_Cols, num.iter, initParam,
                     select_by, fix_Q1, fix_Q2) {
    idx_r <- sort(sample.int(n1, m1, replace = FALSE))
    idx_c <- sort(sample.int(n2, m2, replace = FALSE))
    Xb <- dataMatrix[idx_r, idx_c, drop = FALSE]
    
    local_Q1 <- if (!is.null(fix_Q1)) fix_Q1 else sbmSize$Q1
    local_Q2 <- if (!is.null(fix_Q2)) fix_Q2 else sbmSize$Q2
    local_sbmSize <- list(Q1 = local_Q1, Q2 = local_Q2,
                          explor = if (is.null(sbmSize$explor)) 1.5 else sbmSize$explor)
    
    fit_list <- fitNobiSBM(
      dataMatrix = Xb,
      model = model,
      sbmSize = local_sbmSize,
      exclusive_Rows_Cols = exclusive_Rows_Cols,
      filename = NULL,
      num.iter = num.iter,
      initParam = initParam,
      nbCores = 1L
    )
    best_fit <- .pick_solution(fit_list, select_by = select_by, Q1 = fix_Q1, Q2 = fix_Q2)
    labs <- .get_hard_labels(best_fit)
    
    list(
      idx_r = idx_r,
      idx_c = idx_c,
      labels_row = labs$row,
      labels_col = labs$col,
      Q1 = labs$Q1, Q2 = labs$Q2, ICL = labs$ICL
    )
  }
  
  # Run either with pbapply progress (parallel or serial) or plain lapply/mclapply
  results <- NULL
  if (use_pb) {
    # Serial with progress 
    results <- pbapply::pblapply(
      runs,
      function(b) worker(b, dataMatrix, n1, n2, m1, m2,
                         model, sbmSize, exclusive_Rows_Cols, num.iter, initParam,
                         select_by, fix_Q1, fix_Q2)
    )
  } else {
    # Fallback without pbapply progress
    if (.Platform$OS.type == "windows") {
      results <- lapply(
        runs,
        function(b) worker(b, dataMatrix, n1, n2, m1, m2,
                           model, sbmSize, exclusive_Rows_Cols, num.iter, initParam,
                           select_by, fix_Q1, fix_Q2)
      )
    } else {
      results <- parallel::mclapply(
        runs,
        function(b) worker(b, dataMatrix, n1, n2, m1, m2,
                           model, sbmSize, exclusive_Rows_Cols, num.iter, initParam,
                           select_by, fix_Q1, fix_Q2),
        mc.cores = 1
      )
    }
  }
  
  # Aggregate back on the main process (fast, safe)
  if (verbose && !use_pb) {
    # simple textual progress for the aggregation step
    cat("Aggregating results...\n")
  }
  for (b in seq_along(results)) {
    r <- results[[b]]
    upd <- .update_consensus(cons_same_row, cons_seen_row, r$idx_r, r$labels_row)
    cons_same_row <- upd$cons_same
    cons_seen_row <- upd$cons_seen
    
    upd <- .update_consensus(cons_same_col, cons_seen_col, r$idx_c, r$labels_col)
    cons_same_col <- upd$cons_same
    cons_seen_col <- upd$cons_seen
    
    if (keep_per_run) per_run[[b]] <- r
  }
  
  # Build consensus
  consensus_row <- matrix(0, n1, n1)
  consensus_col <- matrix(0, n2, n2)
  nz_row <- cons_seen_row > 0
  nz_col <- cons_seen_col > 0
  consensus_row[nz_row] <- cons_same_row[nz_row] / cons_seen_row[nz_row]
  consensus_col[nz_col] <- cons_same_col[nz_col] / cons_seen_col[nz_col]
  diag(consensus_row) <- 1
  diag(consensus_col) <- 1
  
  # Full-data fit for cluster-level summaries
  full_fit_list <- fitNobiSBM(
    dataMatrix = dataMatrix,
    model = model,
    sbmSize = list(
      Q1 = if (!is.null(fix_Q1)) fix_Q1 else sbmSize$Q1,
      Q2 = if (!is.null(fix_Q2)) fix_Q2 else sbmSize$Q2,
      explor = if (is.null(sbmSize$explor)) 1.5 else sbmSize$explor
    ),
    exclusive_Rows_Cols = exclusive_Rows_Cols,
    filename = NULL,
    num.iter = num.iter,
    initParam = initParam,
    nbCores = 1
  )
  best_full <- .pick_solution(full_fit_list, select_by = select_by, Q1 = fix_Q1, Q2 = fix_Q2)
  full_labs <- .get_hard_labels(best_full)
  
  row_stab_by_cluster <- .partition_stability(consensus_row, full_labs$row)
  col_stab_by_cluster <- .partition_stability(consensus_col, full_labs$col)
  
  row_vals <- consensus_row[upper.tri(consensus_row, diag = FALSE)]
  col_vals <- consensus_col[upper.tri(consensus_col, diag = FALSE)]
  
  out <- list(
    consensus_row = consensus_row,
    consensus_col = consensus_col,
    seen_row = cons_seen_row,
    seen_col = cons_seen_col,
    stability_summary = list(
      overall = list(
        row_mean = mean(row_vals, na.rm = TRUE),
        row_median = stats::median(row_vals, na.rm = TRUE),
        col_mean = mean(col_vals, na.rm = TRUE),
        col_median = stats::median(col_vals, na.rm = TRUE)
      ),
      by_cluster = list(
        rows = row_stab_by_cluster,
        cols = col_stab_by_cluster
      )
    ),
    full_fit = best_full,
    settings = list(
      B = B, row_frac = row_frac, col_frac = col_frac,
      select_by = select_by, fix_Q1 = fix_Q1, fix_Q2 = fix_Q2,
      model = model, exclusive_Rows_Cols = exclusive_Rows_Cols
    ),
    per_run = per_run
  )
  class(out) <- c("nobiSBM_stability", class(out))
  out
}


# ---- quick plotting helpers -------------------------------------------------
#' Heatmap of Consensus Co-clustering
#'
#' @description
#' Quick visualization of the row or column consensus matrix. Optionally reorders
#' by hierarchical clustering on \code{1 - consensus}.
#'
#' @param stab A \code{"nobiSBM_stability"} object.
#' @param which \code{"row"} or \code{"col"}; which consensus to plot.
#' @param reorder Logical; reorder rows/cols using \code{hclust} on \code{1 - M}.
#' @param method Clustering linkage passed to \code{hclust}; default \code{"complete"}.
#'
#' @return The (possibly reordered) consensus matrix (invisibly).
#'
#' @examples
#' \dontrun{
#' plot_consensus_heatmap(stab, "row")
#' plot_consensus_heatmap(stab, "col", reorder = TRUE, method = "average")
#' }
#' @export
plot_consensus_heatmap <- function(stab, which = c("row","col"),
                                   reorder = TRUE, method = "complete") {
  which <- match.arg(which)
  M <- if (which == "row") stab$consensus_row else stab$consensus_col
  if (reorder) {
    d <- as.dist(1 - M)
    ord <- stats::hclust(d, method = method)$order
    M <- M[ord, ord, drop = FALSE]
  }
  op <- par(no.readonly = TRUE); on.exit(par(op))
  image(t(apply(M, 2, rev)), axes = FALSE, main = paste("Consensus (", which, ")", sep=""))
  box()
  invisible(M)
}

#' Cut a Consensus Matrix into Hard Clusters
#'
#' @description
#' Converts a consensus matrix to a hard partition by hierarchical clustering
#' on \code{1 - consensus} with average linkage.
#'
#' @param stab A \code{"nobiSBM_stability"} object.
#' @param which \code{"row"} or \code{"col"}.
#' @param k Integer; number of clusters to cut.
#'
#' @return Integer vector of cluster labels of length \eqn{n_1} (rows) or \eqn{n_2} (cols).
#'
#' @examples
#' \dontrun{
#' row_part <- cut_consensus(stab, "row", k = stab$full_fit$sbmParam$Q1)
#' }
#' @export
cut_consensus <- function(stab, which = c("row","col"), k) {
  which <- match.arg(which)
  M <- if (which == "row") stab$consensus_row else stab$consensus_col
  d <- stats::as.dist(1 - M)
  stats::cutree(stats::hclust(d, method = "average"), k = k)
}

# obtain bicluster consensus from a stability selection object
get_bicluster_consensus_from_stab <- function(
    stab,
    min_den_row = 1L,
    min_den_col = 1L,
    weight_by_den = TRUE
) {
  C_row <- stab$consensus_row   # n x n
  C_col <- stab$consensus_col   # m x m
  N_row <- stab$seen_row        # same dim as C_row
  N_col <- stab$seen_col        # same dim as C_col
  
  full_fit <- stab$full_fit
  
  # ---------- get biclusters from full fit ----------
  if ("biclusters" %in% names(full_fit)) {
    B_list <- full_fit$biclusters
  } else {
    # fall back: reconstruct from row/col labels, using your existing helper
    labs <- .get_hard_labels(full_fit)  # your helper from A_fun_nobiSBM
    Z1 <- labs$row
    Z2 <- labs$col
    Q1 <- labs$Q1
    Q2 <- labs$Q2
    
    B_list <- list()
    idx <- 0L
    for (q1 in seq_len(Q1)) {
      Rk <- which(Z1 == q1)
      if (!length(Rk)) next
      for (q2 in seq_len(Q2)) {
        Ck <- which(Z2 == q2)
        if (!length(Ck)) next
        idx <- idx + 1L
        B_list[[idx]] <- list(rows = Rk, cols = Ck)
      }
    }
  }
  
  K <- length(B_list)
  if (K == 0L) {
    warning("No biclusters in full_fit.")
    return(NULL)
  }
  
  # helper: average consensus within a subset of indices
  avg_consensus_subset <- function(idx, C, N, min_den, weight_by_den) {
    if (length(idx) <= 1L) {
      return(NA_real_)  # not enough to form pairs
    }
    combs <- t(utils::combn(idx, 2))
    
    c_vals <- C[combs]
    n_vals <- N[combs]
    
    keep <- n_vals >= min_den & !is.na(c_vals)
    if (!any(keep)) {
      return(NA_real_)
    }
    
    if (weight_by_den) {
      return(stats::weighted.mean(c_vals[keep], w = n_vals[keep], na.rm = TRUE))
    } else {
      return(mean(c_vals[keep], na.rm = TRUE))
    }
  }
  
  bic_id       <- integer(K)
  bic_nrow     <- integer(K)
  bic_ncol     <- integer(K)
  bic_row_cons <- numeric(K)
  bic_col_cons <- numeric(K)
  bic_cons     <- numeric(K)
  
  for (k in seq_len(K)) {
    Bk <- B_list[[k]]
    Rk <- Bk$rows
    Ck <- Bk$cols
    
    bic_id[k]   <- k
    bic_nrow[k] <- length(Rk)
    bic_ncol[k] <- length(Ck)
    
    row_con <- avg_consensus_subset(Rk, C_row, N_row,
                                    min_den = min_den_row,
                                    weight_by_den = weight_by_den)
    col_con <- avg_consensus_subset(Ck, C_col, N_col,
                                    min_den = min_den_col,
                                    weight_by_den = weight_by_den)
    
    bic_row_cons[k] <- row_con
    bic_col_cons[k] <- col_con
    
    if (is.na(row_con) || is.na(col_con)) {
      bic_cons[k] <- NA_real_
    } else {
      bic_cons[k] <- row_con * col_con
    }
  }
  
  data.frame(
    bicluster = bic_id,
    n_rows    = bic_nrow,
    n_cols    = bic_ncol,
    row_cons  = bic_row_cons,
    col_cons  = bic_col_cons,
    bic_cons  = bic_cons
  )
}

#' Compute Chance-Level Bicluster Consensus via Random Permutations
#'
#' This function computes the \emph{chance-level} bicluster consensus
#' for each bicluster in a fitted model by generating random biclusters
#' of the same size and evaluating their expected stability.  
#' 
#' The method parallels the random-baseline normalization strategy of 
#' Lee et al. (2010), and provides the quantity
#' \eqn{c^{(k)}_{\mathrm{rand}}} needed for the normalized bicluster 
#' stability measure
#' \deqn{
#'   \tilde c^{(k)} = 
#'     \frac{c^{(k)}_{\mathrm{bic}} - c^{(k)}_{\mathrm{rand}}}
#'          {1 - c^{(k)}_{\mathrm{rand}}}.
#' }
#'
#' Given a bicluster data frame \code{bic_df}, the function repeatedly:
#' \enumerate{
#'   \item randomly permutes all row and column indices,
#'   \item extracts row and column sets of the same size as each bicluster,
#'   \item computes row and column consensus within these random blocks
#'         using \code{stab$consensus_row} and \code{stab$consensus_col},
#'   \item forms a random bicluster consensus as the product of
#'         row and column averages.
#' }
#' The final output is the average random consensus over all \code{T}
#' random permutations for each bicluster.
#'
#' @param bic_df A data frame describing biclusters, typically produced by
#'   \code{get_bicluster_consensus_from_stab()}.  
#'   Must contain the columns:
#'   \itemize{
#'     \item \code{n_rows}: number of rows in bicluster \eqn{k}
#'     \item \code{n_cols}: number of columns in bicluster \eqn{k}
#'   }
#'
#' @param stab A stability object containing:
#'   \itemize{
#'     \item \code{consensus_row}: an \eqn{n \times n} matrix of row 
#'           co-membership probabilities,
#'     \item \code{consensus_col}: an \eqn{m \times m} matrix of column
#'           co-membership probabilities.
#'   }
#'   These matrices are used to evaluate the consensus of randomly 
#'   generated biclusters.
#'
#' @param T Integer; number of random permutations used to estimate the
#'   chance-level consensus. Default is 50.
#'
#' @return A numeric vector of length equal to the number of biclusters.
#'   Each entry is \eqn{c^{(k)}_{\mathrm{rand}}}, the expected bicluster 
#'   consensus under a random assignment of rows and columns.
#'
#' @details
#' The random biclusters preserve only the \emph{sizes} of the row and column
#' sets, not their identities. Row and column indices are permuted jointly
#' across the entire dataset. This ensures that the resulting random
#' biclusters have the same dimensions as those in the original fit but
#' contain no structure, providing an unbiased estimate of chance alignment.
#'
#' @examples
#' \dontrun{
#' bic_df <- get_bicluster_consensus_from_stab(stab)
#' bic_df$c_rand <- compute_random_bic_cons(bic_df, stab, T = 100)
#' }
#'
#' @seealso
#' \code{get_bicluster_consensus_from_stab()},
#' \code{summarize_model_stability()},
#' Lee et al. (2010), \emph{PLoS Computational Biology}.
#'
#' @export
compute_random_bic_cons <- function(bic_df, stab, reps = 50) {
  n <- nrow(stab$consensus_row)
  m <- nrow(stab$consensus_col)
  
  rand_vals <- lapply(seq_len(reps), function(t) {
    perm_r <- sample(seq_len(n), n)
    perm_c <- sample(seq_len(m), m)
    
    sapply(seq_len(nrow(bic_df)), function(k) {
      Rk <- bic_df$n_rows[k]
      Ck <- bic_df$n_cols[k]
      
      rows_rand <- perm_r[seq_len(Rk)]
      cols_rand <- perm_c[seq_len(Ck)]
      
      # compute row and column consensus in the permuted space
      row_con <- mean(stab$consensus_row[rows_rand, rows_rand])
      col_con <- mean(stab$consensus_col[cols_rand, cols_rand])
      
      row_con * col_con
    })
  })
  
  rowMeans(do.call(cbind, rand_vals))
}

