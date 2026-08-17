library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ComplexHeatmap)

# Read in markers from DRSCDB
cho <- read.delim('data/SC_markers_cho.txt')
fu <- read.delim('data/SC_markers_fu.txt')

# Get markers from Seurat object
seuObj <- readRDS('data/yifang_seuratObj.rds')
Idents(seuObj) <- 'cluster'
seuObj <- subset(seuObj, subset = sample != 'Empty_control') # To run on just A/R
 
markers <- FindAllMarkers(seuObj)
markers <- markers %>% 
  filter(grepl('PM', cluster), p_val_adj < 0.05) %>% 
  mutate(cluster = as.character(cluster))

# Filter
cho <- cho %>% filter(p_value < 0.05, cluster_name == 'PM')
fu <- fu %>% filter(p_value < 0.05, grepl('plasmatocyte', cluster_name))

## Bubble heatmap of enriched GO terms per cluster
# Get GO enrichment
bigdf <- data.frame()
clusters <- unique(markers$cluster)

for (i in 1:length(clusters)) {
  clust <- clusters[[i]]
  
  genes <- markers %>% filter(cluster == clust) %>% pull(gene) %>% unique()
  go_res <- enrichGO(genes, OrgDb = org.Dm.eg.db, keyType = 'SYMBOL', ont = 'BP')
  go_res_df <- go_res %>%
    as.data.frame() %>%
    dplyr::select(Description, p.adjust, Count) %>%
    mutate(cluster = clust) %>%
    arrange(p.adjust) %>%
    slice_head(n = 25)
  
  bigdf <- bigdf %>% rbind(go_res_df)
}

# Get enrichment for Cho and Fu as well
genes_cho <- cho %>% pull(gene_symbol) %>% unique()
go_res <- enrichGO(genes_cho, OrgDb = org.Dm.eg.db, keyType = 'SYMBOL', ont = 'BP')
go_res_df <- go_res %>%
  as.data.frame() %>%
  dplyr::select(Description, p.adjust, Count) %>%
  mutate(cluster = 'Cho_PM') %>%
  arrange(p.adjust) %>%
  slice_head(n = 25)

bigdf <- bigdf %>% rbind(go_res_df)

fu <- fu %>%
  mutate(cluster_name = gsub(' plasmatocyte', '_PM', cluster_name),
         cluster_name = paste0('Fu_', cluster_name))
clusters_fu <- unique(fu$cluster_name)

for (i in 1:length(clusters_fu)) {
  clust <- clusters_fu[[i]]
  
  genes <- fu %>% filter(cluster_name == clust) %>% pull(gene_symbol) %>% unique()
  go_res <- enrichGO(genes, OrgDb = org.Dm.eg.db, keyType = 'SYMBOL', ont = 'BP')
  go_res_df <- go_res %>%
    as.data.frame() %>%
    dplyr::select(Description, p.adjust, Count) %>%
    mutate(cluster = clust) %>%
    arrange(p.adjust) %>%
    slice_head(n = 25)
  
  bigdf <- bigdf %>% rbind(go_res_df)
}

# Filter df; top 25 for PM4
go_terms_top_pm4 <- bigdf %>% filter(cluster == 'PM4') %>% arrange(p.adjust) %>% slice_head(n = 25) %>% pull(Description)

all_clusters <- sort(unique(bigdf$cluster))
df <- bigdf %>% 
  filter(Description %in% go_terms_top_pm4) %>%
  mutate(cluster = factor(cluster, levels = all_clusters))

# Plot
ggplot(df, aes(x = cluster, y = Description)) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)),
             shape = 21, color = 'black', stroke = 0.5) +
  scale_fill_gradient(low = 'white', high = '#3393CC') +
  scale_size_continuous(range = c(1, 6)) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'gray80'),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Cell Type', y = 'Enriched GO Term', title = 'Enriched GO Terms for Cluster Markers - Top 25 PM4 Terms',
       fill = '-log10(adj. p-value)', size = 'Count')
ggsave('results/clusterComparison/top25PM4_enrichment_heatmap_ar.png', width = 12, height = 7)

# Filter as top 5 for all
go_terms_top5_all <- bigdf %>% 
  group_by(cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  pull(Description) %>% 
  unique()
df <- bigdf %>% 
  filter(Description %in% go_terms_top5_all) %>%
  mutate(cluster = factor(cluster, levels = all_clusters))

ggplot(df, aes(x = cluster, y = Description)) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)),
             shape = 21, color = 'black', stroke = 0.5) +
  scale_fill_gradient(low = 'white', high = '#3393CC') +
  scale_size_continuous(range = c(1, 6)) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'gray80'),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(x = 'Cell Type', y = 'Enriched GO Term', title = 'Enriched GO Terms for Cluster Markers - Top 5 per Cluster',
       fill = '-log10(adj. p-value)', size = 'Count')
ggsave('results/clusterComparison/top5all_enrichment_heatmap_ar.png', width = 12, height = 9)


### Heatmap of gene x cluster
# Our clusters
min_p <- min(markers$p_val_adj[markers$p_val_adj > 0], na.rm = TRUE) # Set minimum p-value to handle zeroes

# Convert to matrix with values being adjusted p-values
markers_df <- markers %>% 
  mutate(p_val_adj = ifelse(p_val_adj == 0, min_p, p_val_adj),
         log10p = -log10(p_val_adj)) %>%
  dplyr::select(log10p, cluster, gene) %>%
  pivot_wider(names_from = cluster, values_from = log10p) %>%
  dplyr::select(gene, PM1, PM2, PM3, PM4, PM5, PM6, PM7, PM8, PM9)
markers_df[is.na(markers_df)] <- 0

# Cho clusters
min_p_cho <- min(cho$padj[cho$padj > 0], na.rm = TRUE) 

cho_markers_df <- cho %>%
  filter(padj < 0.05) %>%
  mutate(padj = ifelse(padj == 0, min_p_cho, padj),
         log10p = -log10(padj)) %>%
  dplyr::select(log10p, cluster_name, gene_symbol) %>%
  pivot_wider(names_from = cluster_name, values_from = log10p) %>%
  rename(gene = gene_symbol, Cho_PM = PM)

# Fu clusters
min_p_fu <- min(fu$padj[fu$padj > 0], na.rm = TRUE) 

fu_markers_df <- fu %>%
  filter(padj < 0.05) %>%
  mutate(padj = ifelse(padj == 0, min_p_fu, padj),
         log10p = -log10(padj)) %>%
  dplyr::select(log10p, cluster_name, gene_symbol) %>%
  rename(gene = gene_symbol) %>%
  pivot_wider(names_from = cluster_name, values_from = log10p)

# Combine
big_heatmap <- markers_df %>%
  left_join(cho_markers_df, by = 'gene') %>%
  left_join(fu_markers_df, by = 'gene') %>%
  tibble::column_to_rownames('gene') %>%
  select(10:14,1:9) %>%
  as.matrix()

pheatmap(big_heatmap,
         use_raster = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         main = 'Heatmap of Markers per Cluster',
         name = '-log10(adj. p-value)')

  
