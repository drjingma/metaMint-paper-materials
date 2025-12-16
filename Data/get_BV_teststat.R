## obtain test statistics on BV data using one of three methods
## Cai & Liu's approach (cl), asymptotic t test as used in cor.test() and permutation

rm(list=ls())

library(dplyr)
library(ggplot2)
library(cowplot)
library(purrr)
theme_set(theme_cowplot())
library(pheatmap)
library(phyloseq)
library(tibble)
library(tidyr)

OutFolder <- "../Output/biSBM/"
load("../Data/BV_data/BV.rda")
ASV <- t(as(otu_table(BV1), "matrix"))
metabs <- t(dataset$dat)
rownames(metabs) <- dataset$sample_info$sample_ID 

# 1️⃣ Filter sample_info as you already do
dataset$sample_info <- dataset$sample_info %>%
  filter(!is.na(`Nugent Status`)) %>%
  mutate(group = 1 - 1 * (`Nugent Status` %in% "N"))

# 2️⃣ Extract the list of retained sample IDs
keep_ids <- dataset$sample_info$sample_ID

# 3️⃣ Subset ASV to keep only those samples
# Depending on structure of ASV (rows = samples, columns = taxa)
ASV <- ASV[rownames(ASV) %in% keep_ids, , drop = FALSE]
metabs <- metabs[rownames(metabs) %in% keep_ids, , drop=FALSE]

# Optional: verify perfect alignment
stopifnot(all(rownames(ASV) %in% dataset$sample_info$sample_ID))

# 4️⃣ (Optional) Reorder ASV rows to match sample_info exactly
ASV <- ASV[match(dataset$sample_info$sample_ID, rownames(ASV)), , drop = FALSE]
metabs <- metabs[match(dataset$sample_info$sample_ID, rownames(metabs)), , drop = FALSE]
identical(rownames(ASV), dataset$sample_info$sample_ID)
identical(rownames(metabs), dataset$sample_info$sample_ID)

## Test statistics cl were obtained using a constant 0.01 while test statistics pearson were obtained using a constant 0.1.
X <- malabutils::clr.epsilon(ASV,const = 0.1)
Y <- metabs
plot(colMeans(ASV==0))


## ---- one sample Pearson test from Cai & Liu ---- 
#' @param X,Y each being a univariate vector
one_sample_CL_univariate <- function(X,Y){
  n1 <- length(X)
  empcor <- cor(X,Y,method = 'pearson')
  X.scale <- as.numeric(scale(X))
  Y.scale <- as.numeric(scale(Y))
  Theta <- mean((2 * X.scale * Y.scale - empcor * X.scale^2 - empcor * Y.scale^2)^2)
  dataMatrix <- 2*empcor*sqrt(n1/Theta)
  
  return(list(statistic = dataMatrix, estimate = empcor, sd = Theta,
              p.value = 2*pnorm(abs(dataMatrix),lower.tail = FALSE)))
}

two_sample_CL_univariate <- function(x1, y1, x2, y2,
                                     alternative = c("two.sided","greater","less")){
  alternative <- match.arg(alternative)
  # complete cases within each group
  cc1 <- complete.cases(x1, y1); x1 <- x1[cc1]; y1 <- y1[cc1]; n1 <- length(x1)
  cc2 <- complete.cases(x2, y2); x2 <- x2[cc2]; y2 <- y2[cc2]; n2 <- length(x2)
  if (n1 < 4 || n2 < 4) return(c(p_asym = NA_real_, z_asym = NA_real_,
                                 r1 = NA_real_, r2 = NA_real_))
  
  r1 <- one_sample_CL_univariate(x1,y1)
  r2 <- one_sample_CL_univariate(x2,y2)
  
  z <- 2*(r1$estimate - r2$estimate)/sqrt(r1$sd/n1 + r2$sd/n2)
  
  p  <- switch(alternative,
               "two.sided" = 2 * pnorm(-abs(z)),
               "greater"   = 1 - pnorm(z),
               "less"      = pnorm(z))
  
  c(p_asym = p, z_asym = z, r1 = r1$estimate, r2 = r2$estimate)
}
## --------- Asymptotic two-sample Pearson (t-like) ----------
two_sample_pearson_t <- function(x1, y1, x2, y2,
                                 alternative = c("two.sided","greater","less")) {
  alternative <- match.arg(alternative)
  # complete cases within each group
  cc1 <- complete.cases(x1, y1); x1 <- x1[cc1]; y1 <- y1[cc1]; n1 <- length(x1)
  cc2 <- complete.cases(x2, y2); x2 <- x2[cc2]; y2 <- y2[cc2]; n2 <- length(x2)
  if (n1 < 4 || n2 < 4) return(c(p_asym = NA_real_, z_asym = NA_real_,
                                 r1 = NA_real_, r2 = NA_real_))
  r1 <- suppressWarnings(cor(x1, y1, method = "pearson"))
  r2 <- suppressWarnings(cor(x2, y2, method = "pearson"))
  se <- sqrt( (1 - r1^2)^2 / (n1 - 1) + (1 - r2^2)^2 / (n2 - 1) )
  z  <- (r1 - r2) / se
  p  <- switch(alternative,
               "two.sided" = 2 * pnorm(-abs(z)),
               "greater"   = 1 - pnorm(z),
               "less"      = pnorm(z))
  c(p_asym = p, z_asym = z, r1 = r1, r2 = r2)
}

## --------- Permutation two-sample Pearson ----------
two_sample_corr_perm <- function(x1, y1, x2, y2,
                                 B = 5000, seed = NULL,
                                 alternative = c("two.sided","greater","less")) {
  alternative <- match.arg(alternative)
  cc1 <- complete.cases(x1, y1); x1 <- x1[cc1]; y1 <- y1[cc1]; n1 <- length(x1)
  cc2 <- complete.cases(x2, y2); x2 <- x2[cc2]; y2 <- y2[cc2]; n2 <- length(x2)
  if (n1 < 3 || n2 < 3) return(c(p_perm = NA_real_, z_perm = NA_real_))
  r1 <- suppressWarnings(cor(x1, y1, method = "pearson"))
  r2 <- suppressWarnings(cor(x2, y2, method = "pearson"))
  d_obs <- r1 - r2
  x <- c(x1, x2); y <- c(y1, y2)
  g <- c(rep(1L, n1), rep(2L, n2))
  if (!is.null(seed)) set.seed(seed)
  d_null <- replicate(B, {
    gp <- sample(g, length(g), replace = FALSE)
    suppressWarnings(cor(x[gp == 1L], y[gp == 1L], method = "pearson")) -
      suppressWarnings(cor(x[gp == 2L], y[gp == 2L], method = "pearson"))
  })
  p_two <- (sum(abs(d_null) >= abs(d_obs)) + 1) / (B + 1)
  pval  <- switch(alternative,
                  "two.sided" = p_two,
                  "greater"   = (sum(d_null >= d_obs) + 1) / (B + 1),
                  "less"      = (sum(d_null <= d_obs) + 1) / (B + 1))
  z_eq  <- qnorm(1 - p_two/2) * sign(d_obs)
  c(p_perm = pval, z_perm = z_eq)
}

## --------- Wrapper: compute both p-values across taxa for one metabolite ----------
compute_two_sample_pvals_for_metabolite <- function(X, metab, group,
                                                    taxa_subset = NULL,
                                                    B = 5000, seed = 1,
                                                    method = 'pearson',
                                                    alternative = "two.sided",
                                                    eps_rule = "half-min-positive-per-sample",
                                                    eps_const = 1e-6) {
  stopifnot(length(metab) == nrow(X), length(group) == nrow(X))
  if (is.null(taxa_subset)) taxa_subset <- colnames(X)
  if (is.null(taxa_subset) || length(taxa_subset) == 0)
    stop("Provide taxa_subset or ensure colnames(X) exist.")
  group <- droplevels(factor(group))
  if (nlevels(group) != 2) stop("group must have exactly two levels.")
  g1 <- levels(group)[1]; g2 <- levels(group)[2]
  
  # split by group
  idx1 <- which(group == g1)
  idx2 <- which(group == g2)
  
  # iterate taxa
  res <- map_dfr(taxa_subset, function(tx) {
    xj <- X[, tx]
    
    # asymptotic
    a <- two_sample_pearson_t(x1 = xj[idx1], y1 = metab[idx1],
                                       x2 = xj[idx2], y2 = metab[idx2],
                                       alternative = alternative)

    # permutation
    p <- two_sample_corr_perm(x1 = xj[idx1], y1 = metab[idx1],
                              x2 = xj[idx2], y2 = metab[idx2],
                              B = B, seed = seed, alternative = alternative)
    
    # CL univariate
    cl <- two_sample_CL_univariate(x1 = xj[idx1], y1 = metab[idx1],
                                   x2 = xj[idx2], y2 = metab[idx2],
                                   alternative = alternative)
    tibble(
      taxon = tx,
      n1 = length(idx1), n2 = length(idx2),
      r1 = a["r1"] %>% as.numeric(), 
      r2 = a["r2"] %>% as.numeric(),
      diff = (a["r1"] - a["r2"]) %>% as.numeric(),
      p_asym = a["p_asym"] %>% as.numeric(),
      z_asym = a["z_asym"] %>% as.numeric(),
      p_perm = p["p_perm"] %>% as.numeric(),
      z_perm = p["z_perm"] %>% as.numeric(),
      p_cl   = cl['p_asym'] %>% as.numeric(),
      z_cl   = cl['z_asym'] %>% as.numeric()
    )
  })
  
  # add multiple-testing helpers if you like
  res %>%
    mutate(p_asym_bh = p.adjust(p_asym, method = "BH"),
           p_perm_bh = p.adjust(p_perm, method = "BH"),
           p_cl_bh = p.adjust(p_cl, method = "BH"),
           group1 = g1, group2 = g2)
}

## ---------------- Example usage ----------------
# Suppose: X (samples x taxa), metab (vector), group (factor with 2 levels)
# taxa_subset <- taxa_use <- colnames(X)[1]  # or a subset like sample(colnames(X), 500)
group <- dataset$sample_info$group
dataMatrix_CL <- matrix(0, ncol(X), ncol(Y))
dataMatrix_perm <- matrix(0, ncol(X), ncol(Y))
dataMatrix_asym <- dataMatrix_perm
plt_obj <- list()
for (j in 1:ncol(X)){
# for (j in c(1,2,3,35,45,48)){
  cat('running...',j,'..th column\n')
  out <- compute_two_sample_pvals_for_metabolite(X=Y, metab=X[,j], group,
                                                 taxa_subset = colnames(Y),
                                                 B = 5000, seed = 1,
                                                 alternative = "two.sided")
  dataMatrix_CL[j,] <- out$z_cl
  dataMatrix_perm[j,] <- out$z_perm
  dataMatrix_asym[j,] <- out$z_asym
  
  qq_long <- out %>%
    transmute(
      taxon,
      x = z_perm,
      y = z_cl
    )
  plt_obj <- c(plt_obj, list(qq_long))
  
  # p_qq <- ggplot(qq_long, aes(x = x, y = y)) +
  #   geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.7) +
  #   geom_point(alpha = 0.6) +
  #   labs(
  #     title = "Comparison of asymptotic z-statistics vs permutation-based z",
  #     x = expression(z[perm]),
  #     y = expression(z[cl]),
  #   ) +
  #   theme_minimal()
  # print(p_qq)
  
  # save_plot(p_qq, file=paste0('../Data/BV_Analysis/Figures/qqplot_CL_',j,'.png'))
}

saveRDS(list(cl = dataMatrix_CL,
             perm = dataMatrix_perm,
             asym = dataMatrix_asym),
        file = '../Data/BV_Analysis/BV_teststatistics_correlation.rds')

## QQ plot Visualization for selected taxa
names(plt_obj) <- paste0('Taxon',c('1','2','3','35','45','48'))
combined <- imap_dfr(plt_obj, ~ mutate(.x, source = .y)) #%>%
z_long <- as.data.frame(t(dataMatrix_CL[c(1,2,3,35,45,48),])) 
colnames(z_long) <- paste0('Taxon',c('1','2','3','35','45','48'))
z_long <- z_long %>%
  mutate(metabolite = row_number()) %>%        # optional, just an ID
  pivot_longer(
    cols = -metabolite,
    names_to = "taxon",
    values_to = "z"
  )

# Faceted QQ plots vs standard normal
p_qq <- ggplot(z_long, aes(sample = z)) +
  stat_qq(distribution = qnorm, alpha = 0.7) +
  stat_qq_line(distribution = qnorm, linetype = 2, linewidth = 0.5) +
  facet_wrap(~ taxon, scales = "fixed") +
  labs(
    x = "Theoretical quantiles (N(0,1))",
    y = "Observed statistics",
    title = "QQ-plots of taxon-specific association statistics vs N(0,1)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
save_plot(p_qq, file=paste0('../Data/BV_Analysis/Figures/qqplot_CL.png'),
          base_height = 6,base_width = 7)

