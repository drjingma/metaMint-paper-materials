# script to return final selected model after stability analysis on the BV data
rm(list=ls())


library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot())
library(pheatmap)
library(qvalue)
library(phyloseq)
library(dplyr)
library(reshape2)

source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")
source("biSBM/stability_selection.R")

args <- c('Gauss01','0','pearson','4','1')
print(args)
output_dir <- "../Data/BV_Analysis/"
load(paste0(output_dir,'BV.rda'))
metab_info <- dataset$metab_info
taxa_info <- as.data.frame(tax_table(BV1))
taxa_info <- taxa_info %>%
  dplyr::mutate(Genus = ifelse(Genus %in% 'unclassified', paste0(Family,'_', Genus), Genus))
d <- readRDS(paste0(output_dir, 'BV_teststatistics_correlation.rds'))
model <- as.character(args[1])
exclusive_Rows_Cols <- as.integer(args[2])
type <- as.character(args[3])
Q1 <- as.integer(args[4])
Q2 <- as.integer(args[5])
ICL <- (Q1*Q2>0)
alpha.range <- c(0.001,0.005,0.01,seq(0.05, 1, 0.025))

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


sbmSize = switch(
  ICL + 1,  # because switch() in numeric mode is 1-indexed
  list(Q1 = 1:15, Q2 = 1:15, explor = 1.5),       # case ICL == 0
  list(Q1 = Q1, Q2 = Q2, explor = 1.5)            # case ICL == 1
)
nbCores <- 1
modelFamily <- 'Gauss'
n1 <- nrow(dataMatrix)
n2 <- ncol(dataMatrix)

source("biSBM/FDR-R-Code/adaptZ.func.R")
prep <- adaptZ_prepare(as.vector(dataMatrix), gamma = 0.1, model = "Gauss01")
ASC <- lapply(alpha.range, function(a){
  adaptiveZ <- adaptZ_apply(prep,a)
  Ahat <- matrix(0,n1,n2)
  Ahat[adaptiveZ$re] <- 1
  Ahat
}
)
sapply(ASC, sum)/(n1*n2)


res <- fitNobiSBM(dataMatrix, model=model,
                  exclusive_Rows_Cols = exclusive_Rows_Cols,
                  sbmSize = sbmSize,
                  filename = paste0(file.path, '/',file.ending,'.rda'),
                  num.iter = 10,
                  initParam = list(nbOfbeta=1, nbOfPointsPerbeta=NULL,
                                   maxNbOfPasses=3, minNbOfPasses=1),
                  nbCores=parallel::detectCores())
currentSolution <- res[[1]]
gnew <- lapply(alpha.range,function(a)
  graphInferenceNobiSBM(dataMatrix,
                        currentSolution$clustering_row,
                        currentSolution$clustering_col,
                        currentSolution$theta,
                        alpha = a))

Ahat.new <- lapply(gnew, function(g) g$A)
sapply(Ahat.new, sum)/(n1*n2)
cat("ICL...", res[[1]]$sbmParam$ICL, '\n')
## BH procedure ----
ABH <- lapply(alpha.range, function(a) graphInferenceBH(dataMatrix, a)$A)
# sapply(ABH, sum)/(n1*n2)

## qvalue ----
AST <- lapply(alpha.range, function(a) graphInferenceStorey(dataMatrix, a)$A)
# sapply(AST,sum)/(n1*n2)

df <- data.frame(alpha=alpha.range,
                 BH=sapply(ABH, sum)/(n1*n2),
                 Storey=sapply(AST, sum)/(n1*n2),
                 SC=sapply(ASC, sum)/(n1*n2),
                 New=sapply(Ahat.new, sum)/(n1*n2))
plt.power <- reshape2::melt(df,id.vars='alpha') %>%
  dplyr::mutate(Method=variable) %>%
  ggplot(aes(x=alpha,y=value,color=Method,lty=Method)) +
  # scale_shape_manual(values=c(1,0,2,3))+
  # scale_linetype_manual(values = c(2,3,4,1,5)) +
  # geom_point() +
  geom_line() +
  # geom_vline(xintercept = alpha.range,lty=2,color='grey') +
  ylab('Percent of discovered edges') +
  theme(#legend.position = c(0.1, 0.8),
    legend.position = "right",
    axis.title = element_text(size=12),
    legend.text = element_text(size=12),
    legend.title=element_text(size=12)
  )
print(plt.power + xlim(c(0,0.4)) + ylim(c(0,0.5)))

# visualize the density ----
dat <- rbind(data.frame(type='obs',x=as.vector(dataMatrix)),
             data.frame(type='est',x=as.vector(rnbisbm(n1, n2, currentSolution$theta, modelFamily = 'Gauss')$dataMatrix)))

plt.density <- ggplot(dat, aes(x = x)) +
  geom_histogram(data=dat[which(dat$type%in%"obs"),], aes(y = after_stat(density)),  fill = "gray", binwidth = 0.1, alpha=1) +
  stat_function(fun = dnorm, aes(color = 'N(0,1)'),lwd=1) +
  stat_density(data=dat[which(dat$type%in%"est"),], aes(x=x, colour="Fitted"),
               geom="line",position="identity",lty=2,lwd=1) +
  theme(legend.position.inside = c(0.75, 0.75),
        legend.position = "inside")  + labs(color='Density')
plt.density

## Proportion of shared edges ----
names(alpha.range) <- as.character(alpha.range)
prop <- sapply(seq_len(length(alpha.range)), function(a) {
  # Compare overlap with the new method
  # fraction of discoveries in each benchmark
  sapply(list(ABH[[a]], AST[[a]], ASC[[a]]), function(m) {
    overlap <- length(intersect(which(Ahat.new[[a]]==1), which(m==1)))
    overlap / max(length(which(m==1)),1)   
  })     
})
colnames(prop) <- as.character(alpha.range)
rownames(prop) <- c("BH",'Storey','SC')

library(tidyr)
prop_long <- prop %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var='method') %>%
  mutate(method = factor(method, levels = c("BH",'Storey','SC'))) %>%
  pivot_longer(-method, names_to = "FDR", values_to = "prop") %>%
  mutate(FDR = as.numeric(FDR)) %>% 
  filter(FDR <=0.4)

plt_prop <- ggplot(prop_long, aes(x = FDR, y = prop, color = method, lty = method, shape = method)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4")) +
  scale_x_continuous(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)) +
  labs(
    x = "Nominal FDR level",
    y = "Proportion of shared discoveries",
    color = "Method",
    lty = 'Method',
    shape = 'Method'
    # title = "Performance across FDR thresholds"
  ) +
  theme_minimal(base_size = 12) 
# theme(
#   legend.position = "top",
#   panel.grid.minor = element_blank()
# ) 

save_plot(plt_prop,
          file=paste0(file.path,'/', file.ending, "_prop.png"),
          base_asp = 2,base_height = 5)

## define biSBM network ----
id2use <- 6 # corresponds to nominal level 0.1
Ahat.bisbm <- Ahat.new[[id2use]]
rownames(Ahat.bisbm) <- taxa_info$Genus
colnames(Ahat.bisbm) <- metab_info$metabolite
rownames(dataMatrix) <- taxa_info$Genus
colnames(dataMatrix) <- metab_info$metabolite
for (q in 1:max(currentSolution$clustering_row)){
  for (l in 1:max(currentSolution$clustering_col)){
    Ahat.bisbm[currentSolution$clustering_row==q, currentSolution$clustering_col==l] <- 
      currentSolution$theta$nu$mean[q,l] * Ahat.bisbm[currentSolution$clustering_row==q, currentSolution$clustering_col==l]
  }
}
order.row <- order(currentSolution$clustering_row)
order.col <- order(currentSolution$clustering_col)

d_col <- as.dist(1 - cor(dataMatrix))  # correlation distance
hc_col <- hclust(d_col, method = "ward.D2")
order.col <- hc_col$order

mat_ordered <- Ahat.bisbm[order.row, order.col]
datamat_ordered <- dataMatrix[order.row,order.col]
row_clusters <- currentSolution$clustering_row
# col_clusters <- currentSolution$clustering_col

# Step 4: Convert to long format
df <- melt(mat_ordered)
colnames(df) <- c("Genus", "Metabolite", "Value")

df_data <- melt(datamat_ordered)
colnames(df_data) <- c("Genus", "Metabolite", "Value")

df_sub_network <- melt(Ahat.bisbm)
colnames(df_sub_network) <- c("Genus", "Metabolite", "Value")

# Match factor levels to clustering order
df$Genus <- factor(df$Genus, levels = rownames(mat_ordered))
df$Metabolite <- factor(df$Metabolite, levels = colnames(mat_ordered))

# Plot with row and column names in the order estimated
plt_biSBM_network_w_names <- ggplot(df_sub_network, aes(x = Metabolite, y = Genus, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) 

# Plot without row and column names in the estimated order
plt_biSBM_network <- ggplot(df, aes(x = Metabolite, y = Genus, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) 


# Plot without row and column names in the estimated order
plt_data <- ggplot(df_data, aes(x = Metabolite, y = Genus, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) 

A <- mat_ordered
# colnames(A) <- rep("",ncol(A))
# # colnames(A)[setdiff(seq_len(ncol(A)),c(24,70,101,105))] <- NA
# colnames(A)[c(24,70,101,105)] <- colnames(mat_ordered)[c(24,70,101,105)]
A <- A[(rowMeans(abs(A))>0), (colMeans(abs(A))>0)]

# A <- Ahat.bisbm[currentSolution$clustering_row %in% c(2,3),]
df_sub_network <- melt(A)
colnames(df_sub_network) <- c("Genus", "Metabolite", "Value")
# Choose which columns to label
label_cols <- c("GABA",
                "GHB" ,
                "2-hydroxyisovalerate",
                "succinate" ,'tyrosine')

# Convert metabolite to a factor in the desired order
df_sub_network$Metabolite <- factor(df_sub_network$Metabolite, levels = colnames(A))

# Create custom labels: keep only selected, others blank
metabolite_labels <- ifelse(levels(df_sub_network$Metabolite) %in% label_cols,
                            levels(df_sub_network$Metabolite), "")

plt_biSBM_network_w_names <- ggplot(df_sub_network, aes(x = Metabolite, y = Genus, fill = Value)) +
  geom_tile() +
  scale_x_discrete(labels = metabolite_labels) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, limits = c(-0.9, 0.9),
                       oob = scales::squish) +
  labs(x = "Metabolite") +  # make it explicit
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.margin  = margin(t = 5, r = 5, b = 28, l = 5)  # extra bottom space
  )

combined <- plot_grid(
  plot_grid(plt.density, plt.power,
            nrow = 1,
            labels = c("A", "B"),
            scale = 0.95),
  plt_biSBM_network_w_names,
  nrow = 2,
  labels = c("", "C"),
  rel_heights = c(1, 2.5)
)

save_plot(
  filename = paste0(file.path, "/", file.ending, "_power_network.png"),
  plot = combined,
  base_height = 10,
  base_asp = 1
)

## 
save_plot(
  filename = paste0(file.path, "/", file.ending, "_bipartite_network.png"),
  plot = plt_biSBM_network_w_names,
  base_height = 10,
  base_asp = 1
)


