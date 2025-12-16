rm(list=ls())
library(ggplot2)
library(cowplot)
library(dplyr)
library(purrr)
library(RColorBrewer)
theme_set(theme_cowplot())

# ---------- CONFIG ----------
# Provide multiple args as a list of character vectors (each length 15, just like your original `args`)
args_list <- list(
  args <- c("1",'100',"150","200","3","3","1","1","4","8","0.1","0.2",'Gamma','Exp1Gamma', '0'),
  args <- c("1",'100',"150","200","3","3","1","1","4","8","0.1","0.2",'Gamma','Exp1Gamma', '1'),
  args <- c("1",'100',"150","200","3","3","1","0.5","4","8","0.1","0.2",'Gamma','Exp1Gamma', '0'),
  args <- c("1",'100',"150","200","3","3","1","0.5","4","8","0.1","0.2",'Gamma','Exp1Gamma', '1'),
  args <- c("1",'100',"150","200","3","3","1","0.25","1","1","0.1","0.2",'Gamma','Exp1', '0'),
  args <- c("1",'100',"150","200","3","3","1","0.25","1","1","0.1","0.2",'Gamma','Exp1', '1')
  # args <- c("1",'100',"150","200","3","3","1","1","3","6","0.1","0.8",'Gamma','Exp1Gamma', '0'),
  # args <- c("1",'100',"150","200","3","3","1","1","3","6","0.1","0.8",'Gamma','Exp1Gamma', '1'),
  # args <- c("1",'100',"150","200","3","3","1","1","4","8","0.1","0.8",'Gamma','Exp1Gamma', '0'),
  # args <- c("1",'100',"150","200","3","3","1","1","4","8","0.1","0.8",'Gamma','Exp1Gamma', '1'),
  # args <- c("1",'100',"150","200","3","3","1","0.5","4","8","0.1","0.8",'Gamma','Exp1Gamma', '0'),
  # args <- c("1",'100',"150","200","3","3","1","0.5","4","8","0.1","0.8",'Gamma','Exp1Gamma', '1'),
  # args <- c("1",'100',"150","200","3","3","1","0.25","1","1","0.1","0.8",'Gamma','Exp1', '1'),
  # args <- c("1",'100',"150","200","3","3","1","0.25","1","1","0.1","0.8",'Gamma','Exp1', '0')
)
modelFamily <- as.character(args_list[[1]][13])
ICL         <- as.integer(args_list[[1]][15])
shape          <- c(as.numeric(args_list[[1]][9]),as.numeric(args_list[[1]][10]))
prob        <- c(as.numeric(args_list[[1]][11]), as.numeric(args_list[[1]][12]))

OutFolder <- "../Output/newResults_biSBM/"
output_folder_figure <- file.path(OutFolder, "Figures", paste0(prob, collapse = "_"))
if (!dir.exists(output_folder_figure)){
  dir.create(output_folder_figure, showWarnings = TRUE, recursive = TRUE)
}
alpha.range <- c(0.005,0.025,0.05,0.1,0.15,0.25)
myColors <- brewer.pal(9,"Set1")[c(1:5,9)]
colScale <- scale_colour_manual(name = "alpha", values = myColors)

# ---------- HELPERS ----------
safe_read <- purrr::safely(readRDS, otherwise = NULL)

.extract_one <- function(rep_obj) {
  alphas <- setdiff(names(rep_obj), c("Q","randIndex"))
  bind_rows(lapply(alphas, function(a) {
    tab <- rep_obj[[a]]$fdr
    tibble(
      alpha  = as.numeric(a),
      Method = colnames(tab),
      FDR    = as.numeric(tab["FDR", ]),
      TDR    = as.numeric(tab["TDR", ])
    )
  }))
}

.extract_prop <- function(rep_obj) {
  alphas <- setdiff(names(rep_obj), c("Q","randIndex"))
  bind_rows(lapply(alphas, function(a) {
    tibble(
      alpha  = as.numeric(a),
      Method = c("BH","Storey","SC"),
      prop   = rep_obj[[a]]$prop
    )
  }))
}

# --- NEW: extract randIndex (2-element vector) ---
.extract_randIndex <- function(rep_obj) {
  tibble(
    row = rep_obj$randIndex[1],
    col = rep_obj$randIndex[2]
  )
}

.extract_Q <- function(rep_obj) {
  tibble(
    Q = rep_obj$Q
  )
}

# Build a short label for facet strips (customize as you like)
label_from_args <- function(args) {
  # sprintf("prob=%s_%s, ICL=%s, shape=%s_%s, rate=%s",
  sprintf("ICL=%s, shape=%s_%s, rate=%s",
          # args[11], args[12],
          args[15], args[9], args[10], args[8])
}

# Run one setting (one args vector), return summary data frames + label
run_one_setting <- function(args) {
  run_rep <- as.integer(args[1])
  nreps <- as.integer(args[2])
  reps    <- seq_len(nreps)     
  n1 <- as.integer(args[3])
  n2 <- as.integer(args[4])
  Q1 <- as.integer(args[5])
  Q2 <- as.integer(args[6])
  nu0 <- c(1,as.numeric(args[7])) # shape + rate parameter under the null
  nu <- as.numeric(args[8]) # rate parameter under the alternative
  shape_diag <- as.numeric(args[9])
  shape_offdiag <- as.numeric(args[10])
  prob <- c(as.numeric(args[11]),
            as.numeric(args[12]))
  modelFamily <- as.character(args[13])
  model <- as.character(args[14])
  ICL <- as.integer(args[15]) # ICL==2 or 3: indicating misleading the model to incorrect B1,B2
  
  output_folder <- file.path(OutFolder, modelFamily, model, paste0(prob, collapse = "_"))
  file.beginning <- paste0("/res_n1_",n1,"_n2_",n2,"_Q1_",Q1,"_Q2_",Q2, '_shape', 
                           shape_diag, shape_offdiag, 
                           '_rate', nu0[2], nu)
  file.ending <- paste0("_nreps_", nreps, "_ICL", ICL, ".rds")
  
  paths <- sprintf("%s/%s_rep%s%s", output_folder, file.beginning, reps, file.ending)
  res   <- map(paths, safe_read)
  objs  <- compact(map(res, "result"))
  
  # quick report
  fails <- map_chr(res, ~ if (is.null(.x$error)) NA_character_ else conditionMessage(.x$error))
  if (any(!is.na(fails))) {
    message("Failed reads:\n",
            paste(paths[!is.na(fails)], fails[!is.na(fails)], sep=" :: ", collapse="\n"))
  }
  if (!length(objs)) return(NULL)
  
  # combine across replicates
  FDRTDR_long <- imap_dfr(objs, ~ .extract_one(.x) %>% mutate(replicate = as.integer(.y)))
  prop_long   <- imap_dfr(objs, ~ .extract_prop(.x) %>% mutate(replicate = as.integer(.y)))
  rindex      <- imap_dfr(objs, ~ .extract_randIndex(.x) %>% mutate(replicate = as.integer(.y)))
  Qs          <- imap_dfr(objs, ~ .extract_Q(.x) %>% mutate(replicate = as.integer(.y)))
  
  fdr.all.group <- FDRTDR_long %>%
    mutate(alpha=factor(alpha),
           Method=factor(Method, levels=c("New","BH","Storey","SC"))) %>%
    reframe(mfdr = mean(FDR, na.rm = TRUE),
            mtdr = mean(TDR, na.rm = TRUE),
            n    = dplyr::n(),
            .by  = c(Method, alpha))
  
  list(
    label = label_from_args(args),
    fdrtdr = fdr.all.group,
    prop   = prop_long,
    rindex = rindex,
    Q = Qs,
    meta   = list(file.beginning=file.beginning, file.ending=file.ending,
                  modelFamily = modelFamily, model=model,
                  shape = c(shape_diag, shape_offdiag),
                  nu = c(nu, nu0),
                  prob=prob)
  )
}

# ---------- MULTI-SETTINGS FDR plot ----------

# Run all settings, drop NULLs
all_runs <- compact(map(args_list, run_one_setting))

# Combine FDR/TDR summaries with a facet label
fdr_all <- bind_rows(lapply(all_runs, function(z) mutate(z$fdrtdr, Setting = z$label)))

# Plot: one panel per Setting
names(myColors) <- levels(fdr_all$alpha)
vline.df <- tibble(
  alpha = levels(fdr_all$alpha),
  xintercept = as.numeric(as.character(levels(fdr_all$alpha)))
)
# After you build fdr_all:
fdr_all <- fdr_all %>% mutate(alpha = factor(alpha))  # ensure factor

# Dynamic palette sized to the number of alpha levels
n_alpha <- nlevels(fdr_all$alpha)
pal <- if (n_alpha <= 12) {
  RColorBrewer::brewer.pal(max(3, n_alpha), "Set1")[c(1:5,9)]
} else {
  grDevices::hcl.colors(n_alpha, "Dynamic")
}

colScale <- scale_colour_manual(name = "alpha", values = pal)


# Make vlines NOT use the manual color scale
plt_multi <- fdr_all %>%
  dplyr::filter(Method %in% c('New','BH','Storey')) %>%
  # mutate(Setting = recode(Setting,
  #                         "shape=1_3, prob=0.1_0.8, ICL=0" = "B1 = B2 = 3",
  #                         "shape=1_3, prob=0.1_0.8, ICL=1" = "ICL",
  #                         "shape=1_3, prob=0.1_0.8, ICL=2" = 'B1 = B2 = 1',
  #                         "shape=1_3, prob=0.1_0.8, ICL=3" = "B1 = B2 = 5")) %>%
  ggplot(aes(x = mfdr, y = mtdr, shape = Method)) +
  geom_line(alpha = 0.2) +
  geom_point(aes(color = alpha), size = 2.7) +
  scale_shape_manual(values = c(3, 1, 0, 2)) +
  geom_vline(data = vline.df,
             aes(xintercept = xintercept, color = alpha), linetype = "dashed") +
  colScale +
  xlab("FDR") + ylab("TDR") +
  # facet_wrap(~ Setting, scales = "free", ncol = 2,
             # labeller = as_labeller(new_labels)) +
facet_wrap(~ Setting, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 12),
        legend.position = "right")


print(plt_multi)
# Save the multi-panel plot
outfile <- file.path(output_folder_figure, 
                     paste0(modelFamily, '_ICL', ICL, '_shape', paste0(shape, collapse = '_'), 
                            '_',paste0(prob, collapse = '_'),"_fdrtdr_multi_settings.png"))
save_plot(plt_multi, file = outfile, base_height = 5.5, base_width = 11)

# ----------------- MULTI-PANEL PROPORTIONS PLOT -----------------
# Combine 'prop' from each setting and add a facet label
prop_all <- bind_rows(lapply(all_runs, function(z) {
  mutate(z$prop,
         Setting = z$label,
         alpha   = factor(alpha),
         Method  = factor(Method,
                          levels = c("BH","Storey","SC"),
                          labels = c("BH","Storey","SC")))
}))

# If there's nothing to plot, skip gracefully
if (nrow(prop_all) == 0) {
  warning("No proportion data available to plot.")
} else {
  
  # --- Clustered boxplot of proportions ---
  plt_prop_multi <- prop_all %>%
    dplyr::filter(Method %in% c('New','BH','Storey')) %>%
    ggplot(aes(x = alpha, y = prop,
               fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.7) +
    facet_wrap(~ Setting, scales = "free", ncol = 3) +
    scale_fill_brewer(palette = "Set2", name = "Method") +
    xlab(expression(alpha)) + ylab("Proportion of Discoveries") +
    # ggtitle("Proportion of Discoveries by Method and Setting") +
    theme_cowplot() +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # print(plt_prop_multi)
  
  # Save the multi-panel proportions plot
  outfile_prop <- file.path(output_folder_figure, 
                            paste0(modelFamily, '_ICL', ICL, '_shape', paste0(shape, collapse = '_'), 
                                   '_',paste0(prob, collapse = '_'), "_prop_multi_settings.png"))
  save_plot(
    plt_prop_multi,
    file = outfile_prop,
    base_height = 6  # tweak if you have many Settings or alpha levels
    # base_width  = 10
  )
}

# ---- clustering performance ----
rindex <- bind_rows(lapply(all_runs, function(z) {
  mutate(z$rindex,
         Setting = z$label)})) %>%
  group_by(Setting) %>%
  summarise(
    row = mean(row, na.rm = TRUE),
    col = mean(col, na.rm = TRUE),
    # n = n(),
    .groups = "drop"
  )
  
cat('ICL performance in selecting the number of clusters ....\n')
print(rindex)

Qhat <- bind_rows(lapply(all_runs, function(z) {
  mutate(z$Q,
         Setting = z$label)})) %>%
  group_by(Setting) %>%
  summarise(
    mean = mean(Q, na.rm = TRUE),
    sd = sd(Q, na.rm = TRUE),
    # n = n(),
    .groups = "drop"
  )
print(Qhat)

message("Saved: ", outfile)

xtable::xtable(rindex %>% 
                 left_join(Qhat, by = 'Setting')
)
