## results ----
rm(list=ls())
library(ggplot2)
library(cowplot)
library(pheatmap)
theme_set(theme_cowplot())
library(dplyr)
library(purrr)
library(RColorBrewer)
OutFolder <- "../Output/biSBM/"

# biSBM vs SBM varying mean----
model <- "Gauss01"
mu.range <- seq(1,4,0.5)
df <- list()
nreps <- 1:100
for (i in 1:length(mu.range)){
  mu <- mu.range[i]
  file.end <- paste0("20250611_n_100_Q1_2_Q2_3_mu", mu, "_ICL0_bisbm")
  FDRTDR.list = lapply(nreps, function(jid) readRDS(
    paste0(OutFolder,'fdrtdr_',file.end,'_',jid,'.rds')
  ))
  
  df[[i]] <- rbind(data.frame(Method="SBM", replicate = nreps, t(sapply(FDRTDR.list, function(a) c(unlist(a$SBM[,4]), 'time'=as.numeric(a$time[2]))))),
                   data.frame(Method="biSBM", replicate = nreps, t(sapply(FDRTDR.list, function(a) c(unlist(a$biSBM[,4]), 'time' = as.numeric(a$time[1]))))))

}

library(purrr)
df_all <- list_rbind(df, names_to = 'mu') %>%
  mutate(mu = mu.range[mu]) %>%
  group_by(Method, mu) %>% # our group
  summarise(FDR = mean(FDR), TDR = mean(TDR), time = mean(time))
long_df <- df_all %>%
  reshape2::melt(id.vars = c(1,2), 
                 variable.name = "Metric", 
                 value.name = "Value") #%>%
# mutate(Method = factor(Method, levels = c('biSBM','SBM')))

plt <- long_df %>%
  ggplot(aes(y=Value,x=mu,color=Method, shape = Method)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.1, linetype=2) +
  facet_wrap(~ Metric, scales = "free") +
  xlab(expression(mu)) + ylab('Metric') +
  scale_colour_manual(values = c("#377EB8","#E41A1C")) +
  scale_shape_manual(values=c(1,2))

plt

save_plot(plt, file=paste0('../Output//Figures/toy_example_biSBM_average.png'),
          base_height = 3, base_asp = 3)

save_plot(plt, file=paste0('../Output//Figures/', file.end, '_toy_example.png'),
          base_height = 3, base_asp = 3)

## biSBM vs SBM simulation ----
alpha.range <- c(0.005,0.025,0.05,0.1,0.15,0.25)
nreps <- 1:100
# filenames <- '20250517_'
# filenames <- "20250517_n1_150_n2_200_Q1_2_Q2_2_Gauss01_setting4_mu_"
# filenames <- "20250517_n1_150_n2_200_Q1_3_Q2_3_Gauss01_setting2_"
# filenames <- "fdrtdr_20250515_n_101_Q1_2_Q2_1_mu2_ICL_pa_"
filenames <- "fdrtdr_20250527_n_100_Q1_2_Q2_3_mu1_ICL0_sbm_"
FDRTDR.list = lapply(nreps, function(jid) readRDS(
  paste0(OutFolder, filenames, jid,'.rds')
))
Qhat.list <- table(sapply(FDRTDR.list, function(a) a$Q)) # clustering performance is not ideal
fdr.all <-
  rbind(
    do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='biSBM',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$biSBM)))),
    do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='SBM2',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$SBM.2)))),
    # do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='BH',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$BH)))),
    # do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='Storey',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$Storey)))),
    # do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='SC',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$SC07)))),
    do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='SBM',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$SBM))))
  )
rownames(fdr.all) <- NULL

fdr.all <- fdr.all %>%
  mutate(FDR = as.numeric(FDR)) %>%
  mutate(TDR = as.numeric(TDR))
## average FDR and TDR across replicates
fdr.all.group <- fdr.all %>%
  dplyr::reframe(mfdr = mean(FDR,na.rm = TRUE), mtdr = mean(TDR,na.rm = TRUE), .by=c(method,alpha)) %>%
  mutate(alpha=factor(alpha)) %>%
  mutate(method=factor(method))

myColors <- brewer.pal(9,"Set1")[c(1:5,9)]
names(myColors) <- levels(fdr.all.group$alpha)
colScale <- scale_colour_manual(name = "alpha",values = myColors)

plt.fdrtdr.biSBM <- fdr.all.group %>%
  ggplot(aes(x=mfdr, y=mtdr, shape=method)) +
  geom_line(aes(x=mfdr, y=mtdr),alpha=0.2) +
  geom_point(aes(color=alpha),size=3) +
  scale_shape_manual(values=c(3,1,0,2,16,17))+
  geom_vline(xintercept = alpha.range,
             colour = myColors,
             linetype = "dashed") + colScale +
  xlab('FDR') + ylab('TDR') 
plt.fdrtdr.biSBM

# simulation ----
alpha.range <- c(0.005,0.025,0.05,0.1,0.15,0.25)
plt.list <- list()
# names(plt.list) <- c(5,3,4)
Qhat.list <- plt.list
# reps <- setdiff(1:100, c(58,80))
reps <- 1:100
filenames <- rep(NA,5)
filenames[1] <-  '20250507_n1_150_n2_200_Q1_2_Q2_2_Gauss01_setting5_mu2_ICL_res_'
res.all = lapply(reps, function(jid) readRDS(
  paste0(OutFolder, filenames[1], jid,'.rds')
))
A <- readRDS(file=paste0(OutFolder, '20250507_n1_150_n2_200_Q1_2_Q2_2_Gauss01_setting5_network.rds'))
Q <- sapply(res.all, function (a) unlist(getBestQ(a, F)))
filenames[2] <- '20250507_n1_150_n2_200_Q1_3_Q2_3_Gauss01_setting2_'
filenames[4] <- '20250507_n1_150_n2_200_Q1_2_Q2_2_Gauss01_setting4_mu2_ICL_'
filenames[3] <- '20250507_n1_150_n2_200_Q1_2_Q2_1_Gauss01_setting5_mu2_'
filenames[5] <- '20250507_n1_150_n2_200_Q1_2_Q2_2_Gauss01_setting5_mu2_ICL_'
filenames[4] <- '20250517_'
for (setting in c(2,4,5)){
  setting <- 4
  
  # filenames <- list.files(OutFolder)
  # index <- grep('20250505_n1_150_n2_200_Q1_3_Q2_3_Gauss01_setting3', filenames)
  # # Extract numbers before .rds
  # numbers <- sub(".*_([0-9]+)\\.rds$", "\\1", filenames[index])
  FDRTDR.list = lapply(reps, function(jid) readRDS(
    paste0(OutFolder, filenames[setting], jid,'.rds')
  ))
  Qhat.list[[setting]] <- table(sapply(FDRTDR.list, function(a) a$Q)) # clustering performance is not ideal
  fdr.all <-
    rbind(
      do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='biSBM',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$biSBM)))),
      # do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='SBM',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$SBM)))),
      # do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='SBM2',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$SBM.2)))),
      do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='BH',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$BH)))),
      do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='Storey',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$Storey)))),
      do.call(rbind,lapply(1:length(FDRTDR.list),function(i) data.frame(method='SC',replicate=i,alpha=alpha.range,t(FDRTDR.list[[i]]$SC07))))
      # do.call(rbind,lapply(1:length(index),function(i) data.frame(method='SC09',replicate=i,alpha=alpha.range,t(FDRTDR[[i]]$SC09))))
    )
  rownames(fdr.all) <- NULL
  
  fdr.all <- fdr.all %>%
    mutate(FDR = as.numeric(FDR)) %>%
    mutate(TDR = as.numeric(TDR))
  ## average FDR and TDR across replicates
  fdr.all.group <- fdr.all %>%
    dplyr::reframe(mfdr = mean(FDR,na.rm = TRUE), mtdr = mean(TDR,na.rm = TRUE), .by=c(method,alpha)) %>%
    mutate(alpha=factor(alpha)) %>%
    mutate(method=factor(method,levels=c('biSBM','BH','Storey','SC'), labels=c('New','BH','Storey','SC')))
  
  myColors <- brewer.pal(9,"Set1")[c(1:5,9)]
  names(myColors) <- levels(fdr.all.group$alpha)
  colScale <- scale_colour_manual(name = "alpha",values = myColors)
  
  plt.fdrtdr.biSBM <- fdr.all.group %>%
    ggplot(aes(x=mfdr, y=mtdr, shape=method)) +
    geom_line(aes(x=mfdr, y=mtdr),alpha=0.2) +
    geom_point(aes(color=alpha),size=3) +
    scale_shape_manual(values=c(3,1,0,2))+
    geom_vline(xintercept = alpha.range,
               colour = myColors,
               linetype = "dashed") + colScale +
    xlab('FDR') + ylab('TDR') 
  plt.list[[setting]] <- plt.fdrtdr.biSBM
}

legend = cowplot::get_plot_component(plt.list[[5]], 'guide-box-right', return_all = TRUE)

plt.list[[5]] <- plt.list[[5]] + theme(legend.position="none")
plt.list[[4]] <- plt.list[[4]] + theme(legend.position="none")
plt.list[[2]] <- plt.list[[2]] + theme(legend.position="none")

plt <- plot_grid(plot_grid(plt.list[[2]], plt.list[[4]], plt.list[[5]], nrow=1,
                           labels = c('a','b','c'), scale = 0.95),
                 legend,
                 nrow = 1,
                 rel_widths = c(0.875,.125))

save_plot(plt, file=paste0(OutFolder,'/Figures/simulation_biSBM.png'),
          base_height = 3.5, base_asp = 3.5)

