library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Read in
mat <- read.delim('data/matrix_GEO.txt', check.names = FALSE, row.names = 1)
meta <- read.delim('data/metadata_GEO.txt')

# Make Seurat
meta <- meta %>% filter(!barcode == '')
rownames(meta) <- meta$barcode
mat <- as.matrix(mat)
all(colnames(mat) == rownames(meta))
seuObj <- CreateSeuratObject(mat, meta.data = meta)

Idents(seuObj) <- 'cluster'

meta <- seuObj@meta.data
meta <- meta %>% mutate(sample = ifelse(sample == 'Empty_control', 'Control', sample))
seuObj@meta.data <- meta

seuObj <- NormalizeData(seuObj)

# Violin plots - one plot for each sample, separated by cell type
samps <- c('Control', 'AML1-ETO', 'RasV12')
plots <- list()

for (i in 1:length(samps)) {
  samp <- samps[[i]]

  seuObj_subset <- subset(seuObj, sample == samp) 
  
  p <- VlnPlot(seuObj_subset, features = 'Mrtf', group.by = 'cluster') + 
    labs(title = paste0(samp, ' Sample - Mrtf Violin Plot'))
  
  plots[[i]] <- p
}

combined_plot <- plots[[1]] / plots[[2]] / plots[[3]]

ggsave('results/violinPlots/Mrtf_by_sample.png', combined_plot, width = 10, height = 12)

# Opposite - one plot for each cell type, separated by sample
clusts <- c('CC1', 'CC2', 'LM1', 'LM2', 'PM1', 'PM2', 
            'PM3', 'PM4', 'PM5', 'PM6', 'PM7') # PM8 and PM9 only present in control; don't include
plots <- list()

for (i in 1:length(clusts)) {
  clust <- clusts[[i]]

  seuObj_subset <- subset(seuObj, cluster == clust) 
  Idents(seuObj_subset) <- 'sample'
  Idents(seuObj_subset) <- factor(Idents(seuObj_subset), levels = c('Control', 'AML1-ETO', 'RasV12'))
  
  p <- VlnPlot(seuObj_subset, features = 'Mrtf') + 
    labs(title = paste0(clust, ' Cells - Mrtf Violin Plot'))
  
  plots[[i]] <- p
}

combined_plot <- wrap_plots(plots, ncol = 4)
ggsave('results/violinPlots/Mrtf_by_cluster.png', combined_plot, width = 18, height = 12)

plots_of_int <- wrap_plots(list(plots[[4]], plots[[7]]), ncol = 2)
ggsave('results/violinPlots/Mrtf_PM3_LM2.png', plots_of_int, width = 12, height = 6)
