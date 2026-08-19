library(monocle3)
library(dplyr)
library(patchwork)
library(ggplot2)

# Read in
mat <- read.delim('data/matrix_GEO.txt', check.names = FALSE, row.names = 1)
meta <- read.delim('data/metadata_GEO.txt')
gene_ann <- read.csv('data/gene_annotations.csv', row.names = 1)

## SPLIT INTO GENOTYPE
samp <- 'AML1-ETO'  # Empty_control, RasV12, AML1-ETO

meta <- meta %>% filter(!barcode == '')
rownames(meta) <- meta$barcode

cells <- meta %>% filter(sample == samp) %>% pull(barcode)

mat <- mat %>% select(all_of(cells))
meta <- meta %>% filter(barcode %in% cells)

# Fix formatting
mat <- as.matrix(mat)

all(colnames(mat) == rownames(meta))
all(rownames(mat) == rownames(gene_ann)) 

# Create object
cds <- new_cell_data_set(mat, cell_metadata = meta, gene_metadata = gene_ann)
cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds, preprocess_method = 'PCA')
colData(cds)@listData[["Celltype"]] <- meta$cluster
cds@clusters@listData[["UMAP"]][["clusters"]] <- meta$cluster
cds@clusters@listData[["UMAP"]][["partitions"]] <- meta$cluster

# Reset UMAP and PCA data
seuObj <- readRDS('data/yifang_seuratObj.rds')
seuObj <- subset(seuObj, subset = sample == samp)
cds@int_colData@listData[["reducedDims"]]@listData[["PCA"]] <- seuObj@reductions[["pca"]]@cell.embeddings
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seuObj@reductions[["umap"]]@cell.embeddings

# run monocle3
cds <- cluster_cells(cds)
cds <- learn_graph(cds,use_partition = F)

plot_cells(cds,
           color_cells_by = "Celltype",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=2.5,
           group_label_size = 3)

# Helper function from monocle3 vignette that chooses the root node as the highest proportion of PM1
get_earliest_principal_node <- function(cds, cluster="PM1"){
  cell_ids <- which(colData(cds)[, "Celltype"] == cluster)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

png(paste0('results/trajectoryAnalysis_v2/', samp, '_trajectory.png'), width = 12, height = 6, units = 'in', res = 300)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, 
           label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=2.5) + 
  plot_cells(cds, color_cells_by = "Celltype", label_cell_groups=TRUE, group_cells_by = 'Celltype',
           label_leaves=FALSE, label_branch_points=FALSE, group_label_size=3.5)
dev.off()

# Get pseudotime table, UMAP coords, and cluster membership to plot with more control
pseudotime_table <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime_table <- pseudotime_table[rownames(meta)] %>% data.frame()
colnames(pseudotime_table) <- 'Pseudotime'
pseudotime_table <- pseudotime_table %>% tibble::rownames_to_column('barcode')

umap_df <- seuObj@reductions[["umap"]]@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column('barcode')

cluster_mem <- meta %>% select(-sample)

# Merge all together
df <- umap_df %>% left_join(cluster_mem, by = 'barcode') %>% left_join(pseudotime_table, by = 'barcode')

pal = c('PM1' = "#7ED321", 'PM2' = "#F5A623", 'PM3' = "#417505", 
           'PM4' = "#F1DF05", 'PM5' = "#50E3C2", 'PM6' = "#4A4A4A", 
           'PM7' = "#4A90E2", 'PM8' = "#9a9a9a", 'PM9' = "#000000",
           'CC1' = "#F199A3", 'CC2' = "#D0021B", 'LM1' = "#BA98FF", 'LM2' = "#BD10E0") 

# Plot
p_time <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = Pseudotime)) +
  geom_point(size = 1) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = 'UMAP by Pseudotime') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
p_clust <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = cluster)) +
  geom_point(size = 1) +
  scale_color_manual(values = pal) +
  theme_minimal() +
  labs(title = 'UMAP by Cell Type', color = 'Cell Type') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

p_time + p_clust

ggsave(paste0('results/trajectoryAnalysis_v2/', samp, '_pseudotime_oldcolors.png'), width = 12, height = 6)
