library(monocle3)
library(dplyr)
library(patchwork)
setwd("~/Documents/scRNA/with/Sudhir/leukemia_reanalysis")

# Read in
mat <- read.delim('data/matrix_GEO.txt', check.names = FALSE, row.names = 1)
meta <- read.delim('data/metadata_GEO.txt')
gene_ann <- read.csv('data/gene_annotations.csv', row.names = 1)

## SPLIT INTO GENOTYPE
meta <- meta %>% filter(!barcode == '')
rownames(meta) <- meta$barcode

cells <- meta %>% filter(sample == 'Empty_control') %>% pull(barcode) # Empty_control, RasV12, AML1-ETO

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
#cds@int_colData@listData[["reducedDims"]]@listData[["PCA"]] <- 
#  PC_DB_DT@reductions[["pca"]]@cell.embeddings     # To use old PCA and UMAP, which I don't have
#cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- PC_DB_DT@reductions[["umap"]]@cell.embeddings
colData(cds)@listData[["Celltype"]] <- meta$cluster
cds@clusters@listData[["UMAP"]][["clusters"]] <- meta$cluster
cds@clusters@listData[["UMAP"]][["partitions"]] <- meta$cluster

# Check UMAP
#genes_check <- c('hop', 'EcR', 'foxo', 'CalpB', 'gal4', 'Mrtf')
genes_check <- c('NimC1', 'eater', 'lz', 'PPO1', 'atilla', 'PPO3')

plot_cells(cds, genes = genes_check, label_groups_by_cluster=FALSE,  color_cells_by = "Celltype", show_trajectory_graph = FALSE)

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

png('results/trajectoryAnalysis/emptyCtrl_trajectory.png', width = 12, height = 6, units = 'in', res = 300)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, 
           label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=2.5) + 
  plot_cells(cds, color_cells_by = "Celltype", label_cell_groups=TRUE, group_cells_by = 'Celltype',
           label_leaves=FALSE, label_branch_points=FALSE, group_label_size=3.5)
dev.off()

png('results/trajectoryAnalysis/emptyCtrl_trajectory_WITHLEGEND.png', width = 12, height = 6, units = 'in', res = 300)
plot_cells(cds, color_cells_by = "Celltype", label_cell_groups=FALSE, group_cells_by = 'Celltype',
             label_leaves=FALSE, label_branch_points=FALSE, group_label_size=3.5)
dev.off()

# Get table
pseudotime_table <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime_table <- pseudotime_table[rownames(meta)] %>% data.frame()


# # SUBSET for specific branch to end nodes
# cds_sub <- choose_graph_segments(cds)
# cds_sub <- preprocess_cds(cds_sub, method = 'PCA')
# cds_sub <- reduce_dimension(cds_sub, preprocess_method = 'PCA')
# cds_sub <- cluster_cells(cds_sub)
# cds_sub <- learn_graph(cds_sub,use_partition = F)
# cds_sub <- order_cells(cds_sub, root_pr_nodes=get_earliest_principal_node(cds_sub))
# plot_cells(cds_sub,
#            color_cells_by = "pseudotime",
#            label_cell_groups=FALSE,
#            label_leaves=TRUE,
#            label_branch_points=FALSE,
#            graph_label_size=2.5)
# plot_cells(cds_sub,
#            color_cells_by = "seurat_clusters",
#            label_cell_groups=FALSE,
#            label_leaves=TRUE,
#            label_branch_points=FALSE,
#            graph_label_size=2.5)

# Below is just copied from trachea paper, maybe look into it later
# identify genes differentially regulated in pseudotime
pr_test_res <- graph_test(cds, neighbor_graph="principal_graph")
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=0.001)
names(gene_module_df)[1] <- "gene_short_name"
pr_merge <- merge(pr_test_res, gene_module_df, by = "gene_short_name")






