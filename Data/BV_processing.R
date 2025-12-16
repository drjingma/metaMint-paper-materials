
rm(list=ls())
library(RColorBrewer)
# display.brewer.all()
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot())
library(pheatmap)
library(qvalue)
library(tidyverse)

source("biSBM/dataGeneration_biSBM.R")
source("biSBM/auxiliaryFunctions_biSBM.R")
source("biSBM/VEMinitialization_biSBM.R")
source("biSBM/VEMalgorithm_biSBM.R")
source("biSBM/testProcedure_biSBM.R")
source("arXiv/lib/censoredGGM.R")
source("arXiv/lib/preprocess_lib.r")


out_dir <- "../Data/BV_data/"

# read data and organize into a phyloseq object
taxa <- read_csv('../Data/BV_data/GregorReid/srep14174-s4-taxa.csv',skip=1)
sampleInfo <- read_csv('../Data/BV_data/GregorReid/srep14174-s1-metadata.csv', skip=1)
otumat <- as.matrix(taxa[,-1])
rownames(otumat) <- paste0('spe',seq(nrow(otumat)))
tnames <- taxa$Taxa
tmp <- lapply(tnames, function(s) unlist(strsplit(s, ";")))
taxmat <- do.call(rbind, tmp)
taxmat <- cbind(taxmat,taxa$Taxa)
colnames(taxmat) <- c('Phylum', 'Class', 'Order', 'Family', 'Genus',"Taxa")
rownames(taxmat) <- rownames(otumat)
taxmat <- taxmat[-49,]

TAX <- tax_table(taxmat)
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
OTU[30,] <- OTU[30,] + OTU[49,]
OTU <- OTU[-49,]

SAM <- as.data.frame(sampleInfo[1:131,])
rownames(SAM) <- SAM$Patient_ID
SAM <- SAM[match(colnames(otumat),rownames(SAM)),]
BV <- phyloseq(OTU, TAX, sample_data(SAM))

prevdf = apply(X = otu_table(BV),
               MARGIN = ifelse(taxa_are_rows(BV), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(BV),
                    tax_table(BV))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# set prevalence threshold
prevalenceThreshold = 0.2 * nsamples(BV)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
BV1 = prune_taxa(keepTaxa, BV)
BV1

# obtain abundance and apply mCLR transformation 
count_OTU <- as(otu_table(BV1), "matrix")
abundance_mclr <- clr_epsilon(t(count_OTU))

# metabolomic data
metab <- read_csv("../Data/BV_data/GregorReid/srep14174-s2-metabolites.csv", skip=2)
metab2 <- read_csv("../Data/BV_data/GregorReid/srep14174-s5-metabolites.csv", skip=2)

dataset <- list()
dataset$sample_info <- as.data.frame(t(metab[1:3,-c(2:7)]), stringsAsFactors = F)
colnames(dataset$sample_info) <- sapply(dataset$sample_info[1,], as.character)
dataset$sample_info <- dataset$sample_info[-1,]
dataset$sample_info$sample_ID <- rownames(dataset$sample_info)
rownames(dataset$sample_info) <- NULL
dataset$sample_info$diversity <- as.numeric(dataset$sample_info$diversity)
dataset$sample_info$Pregnant <- as.numeric(dataset$sample_info$Pregnant)

dataset$metab_info <- data.frame('metab_ID'=seq(1,nrow(metab)-3),'metabolite'=metab$sample_ID[-c(1:3)], stringsAsFactors = FALSE)
rownames(dataset$metab_info) <- NULL

dataset$dat <- as.data.frame(metab[-c(1:3),-c(1:7)])
dataset$dat <- apply(dataset$dat, 2, as.numeric)
colnames(dataset$dat) <- NULL

SAM <- as.data.frame(sample_data(BV1))
dataset <- order.samples(dataset, match(SAM$Patient_ID, dataset$sample_info$sample_ID))
dataset <- order.metabs(dataset, match(sort(dataset$metab_info$metabolite), dataset$metab_info$metabolite))
rownames(dataset$dat) <- paste0('m',1:nrow(dataset$dat))
data.input.MTB <- list(taxa = t(abundance_mclr), metabolites = dataset$dat)

save(dataset,BV1,file=paste0(out_dir,"BV.rda"))

