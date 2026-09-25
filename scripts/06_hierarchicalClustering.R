library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(ComplexHeatmap)

# Read in markers and background gene set (from seurat object)
cho <- read_excel('data/cho_markers_subcluster.xlsx')
tattikota <- read.csv('data/tattikota_2020_markers_all.csv')
markers <- read.csv('data/markers_all.csv')

background_genes <- read.delim('data/backgroundGenes.txt', header = FALSE) %>% pull(V1)

# Get top 100 for each dataset's clusters
cho <- cho %>% filter(grepl("PM|LM|CC", Celltype)) %>% mutate(Celltype = gsub(' ', '', Celltype)) %>% rename(cluster = Celltype, gene = Gene)
top100cho <- cho %>% 
  group_by(cluster) %>%
  slice_max(order_by = avg_l2fc, n = 100) %>%
  ungroup() %>% 
  mutate(cluster = paste0(cluster, '_Cho')) %>%
  select(cluster, gene)

top100tattikota <- tattikota %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100) %>%
  ungroup() %>% 
  mutate(cluster = paste0(cluster, '_Tattikota')) %>%
  select(cluster, gene)

top100markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100) %>%
  ungroup() %>%
  select(cluster, gene)

# Get all markers as a list for each cluster
top100 <- rbind(top100cho, top100tattikota, top100markers) %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = 'drop')
clusters <- unique(top100$cluster)

# Use hypergeometric test to get a p-value of similarity between each cluster
p_vals <- expand.grid(cluster1 = clusters, cluster2 = clusters) %>%
  rowwise() %>%
  mutate(genes1 = list(top100$genes[[match(cluster1, top100$cluster)]]),
         genes2 = list(top100$genes[[match(cluster2, top100$cluster)]]),
         overlap = length(intersect(genes1, genes2)),
         N = length(background_genes),
         m = length(genes1),
         k = length(genes2),
         p_value = phyper(overlap - 1, m, (N - m), k, lower.tail = FALSE),
         similarity = -log10(p_value)) %>%
  ungroup() %>%
  select(cluster1, cluster2, overlap, p_value, similarity)

# Values are too high for visualization; cap similarity score at 100 (only for 50 and 100, not 10)
p_vals <- p_vals %>%
  mutate(similarity = ifelse(similarity > 100, 100, similarity))

# Make a matrix for a heatmap
mat <- p_vals %>%
  select(cluster1, cluster2, similarity) %>%
  pivot_wider(names_from = cluster2, values_from = similarity) %>%
  tibble::column_to_rownames('cluster1') %>%
  as.matrix()

# For hierarchical clustering, need to convert it to distance
max_similarity <- max(mat, na.rm = TRUE)
dist_mat <- max_similarity - mat

# Add annotation bars
ann_df <- data.frame(cluster = clusters) 
ann_df <- ann_df %>%
  mutate(Celltype = ifelse(grepl('CC', cluster), 'Crystal cell', 
                           ifelse(grepl('LM', cluster), 'Lamellocyte', 'Plasmatocyte'))) %>%
  tibble::column_to_rownames('cluster')

annotation_colors <- list(Celltype = RColorBrewer::brewer.pal(n = length(unique(ann_df$Celltype)), name = "Set2"))
names(annotation_colors$Celltype) <- unique(ann_df$Celltype)

png('results/clusterComparison_v2/hierarchicalClustering_10markers.png', width = 10, height = 8, units = 'in', res = 300)
pheatmap(mat,
         annotation_col = ann_df,
         annotation_row = ann_df,
         annotation_colors = annotation_colors,
         clustering_distance_rows = as.dist(dist_mat),
         clustering_distance_cols = as.dist(dist_mat),
         clustering_method = "complete",
         color = colorRampPalette(c('white', 'tomato', 'red'))(99),
         main = "Hierarchical Clustering of Clusters by Top 100 Marker Similarity",
         name = '-log10(P)')
dev.off()

