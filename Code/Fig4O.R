WT_Cellbin_clust <- load("/location_where_you_input_file_is/D80_WTCellbin_semi_out1.RData")
WT_Clust <- semisup
Can_Cellbin_clust <- load("/location_where_you_input_file_is/D80_CanCellbin_semi_out1.RData")
Can_Clust <- semisup

WT <- readRDS("/group/aawanggrp/siyulin/InSitu_Merging/RDS/D80_Half_WT.rds")
Can <- readRDS("/group/aawanggrp/siyulin/InSitu_Merging/RDS/D80_Half_Can.rds")

desired_clustersWT <- "clst5283"
desired_clustersCan <- "clst5283"

# Subset both wt and canavan based on your desired cluster identities
subset_WT_Clust <- WT_Clust$clust[WT_Clust$clust %in% desired_clustersWT]
subset_Can_Clust <- Can_Clust$clust[Can_Clust$clust %in% desired_clustersCan]

# extract the cell ids from cluster files
subset_WT_cell_ids <- names(subset_WT_Clust)
subset_Can_cell_ids <- names(subset_Can_Clust)

# extract the  cell ids from the RDS files
WT_cellbin_cell_ids <- rownames(WT@meta.data)
Can_cellbin_cell_ids <- rownames(Can@meta.data)

# find the matched cell ids from both files
WT_matching_cell_ids <- intersect(subset_WT_cell_ids, WT_cellbin_cell_ids)
WT_number_of_matching_cells <- length(WT_matching_cell_ids)
print(WT_number_of_matching_cells) 

Can_matching_cell_ids <- intersect(subset_Can_cell_ids, Can_cellbin_cell_ids)
Can_number_of_matching_cells <- length(Can_matching_cell_ids)
print(Can_number_of_matching_cells) 

WT_cluster_information <- subset_WT_Clust
Can_cluster_information <- subset_Can_Clust

# store only the information about the matched clusters for further use
WT_matched_clusters <- WT_cluster_information[WT_matching_cell_ids]
Can_matched_clusters <- Can_cluster_information[Can_matching_cell_ids]

# create a new column, name it as "matched_clusters" and fill it up with NA
WT$matched_clusters <- NA 
Can$matched_clusters <- NA 

# fill in the matched cluster information
WT$matched_clusters[WT_matching_cell_ids] <- WT_matched_clusters
head(WT@meta.data)

Can$matched_clusters[Can_matching_cell_ids] <- Can_matched_clusters
head(Can@meta.data)

metadata <- Can@meta.data

# Identify rows (cells) where 'matched_clusters' is not NA
valid_rows <- !is.na(metadata$matched_clusters)

# Subset the Seurat object based on these valid rows
Can_subset <- Can[, valid_rows]

metadata <- WT@meta.data

# Identify rows (cells) where 'matched_clusters' is not NA
valid_rows <- !is.na(metadata$matched_clusters)

# Subset the Seurat object based on these valid rows
WT_subset <- WT[, valid_rows]

WT_subset$orig.ident <- "D80_WT"
Can_subset$orig.ident <- "D80_Can"

D80 <- merge(WT_subset, y = Can_subset)

D80 <- NormalizeData(D80, normalization.method = "LogNormalize", scale.factor = 10000)

D80 <- FindVariableFeatures(D80, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(D80)
D80 <- ScaleData(D80, features = all.genes)

D80 <- RunPCA(D80, features = VariableFeatures(object = D80))
print(D80[["pca"]], dims = 1:5, nfeatures = 5)

de.results_D80_2 <- FindMarkers(D80,
                                ident.1 = "D80_Can", 
                                ident.2 = "D80_WT",
                                group.by = "orig.ident",
                                test.use = "MAST",
                                min.pct = 0.4,
                                logfc.threshold = 0.25)

write.xlsx(de.results_D80_2, file = "D80_5283_04.xlsx", rowNames = TRUE)

genes_of_interest <- c("Pikfyve", "Jmjd6", "Cops5", "Kcnb2", "Srsf2", "Dcaf6", "Lonrf2", "Xkr4", "Rp1", "Scarna6", "Rab17", "Rpe", "Tyw5", "Cflar", "Prdm14", "Ankrd39", "Nefh", "Wdr12", "Tfap2b")

scaled_data <- GetAssayData(D80, slot = "scale.data")[genes_of_interest, ]
D80$Condition <- ifelse(D80$orig.ident == "D80_Can", "D80_Can", "D80_WT")
cells_can <- colnames(scaled_data)[D80$orig.ident == "D80_Can"]
cells_wt <- colnames(scaled_data)[D80$orig.ident == "D80_WT"]
set.seed(200) 
sample_cells_can <- sample(cells_can, size = round(length(cells_can) * 1))
sample_cells_wt <- sample(cells_wt, size = round(length(cells_wt) * 1))
sample_cells <- c(sample_cells_can, sample_cells_wt)
conditions <- c(rep("D80_Can", length(sample_cells_can)), rep("D80_WT", length(sample_cells_wt)))
names(conditions) <- sample_cells
sample_cells_sorted <- sample_cells[order(conditions)]
plotdat <- as.matrix(scaled_data[, sample_cells_sorted])

outlier_cells <- colnames(plotdat)[apply(plotdat, 2, function(cell_data) any(cell_data > 100))]
plotdat <- plotdat[, !(colnames(plotdat) %in% outlier_cells)]
conditions <- conditions[!(names(conditions) %in% outlier_cells)]

col_ha <- HeatmapAnnotation(df = data.frame(Condition = conditions[sample_cells_sorted]))
can_indices_in_plotdat <- which(D80$orig.ident[colnames(plotdat)] == "D80_Can")
wt_indices_in_plotdat <- which(D80$orig.ident[colnames(plotdat)] == "D80_WT")
ordered_indices_in_plotdat <- c(can_indices_in_plotdat, wt_indices_in_plotdat)

gene_order <- c("Pikfyve", "Jmjd6", "Cops5", "Kcnb2", "Srsf2", "Dcaf6", "Lonrf2", "Xkr4", "Rp1", "Scarna6", "Rab17", "Rpe", "Tyw5", "Cflar", "Prdm14", "Ankrd39", "Nefh", "Wdr12", "Tfap2b")

ordered_rows <- match(gene_order, rownames(plotdat))

plotdat <- as.matrix(scaled_data[, sample_cells_sorted])

outlier_cells <- colnames(plotdat)[apply(plotdat, 2, function(cell_data) any(cell_data > 100))]
plotdat <- plotdat[, !(colnames(plotdat) %in% outlier_cells)]

conditions <- conditions[!(names(conditions) %in% outlier_cells)]
sample_cells_sorted <- sample_cells_sorted[!(sample_cells_sorted %in% outlier_cells)]

col_ha <- HeatmapAnnotation(df = data.frame(Condition = conditions[sample_cells_sorted]))

breakpoints <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)
custom_colors <- c("#001260FF", "#014683FF", "#327FA7FF", "#CCDFE9FF", "#E9ECEBFF", "#EAE0CDFF", "#DAC9A5FF", "#CAB17EFF", "#CB8E62", "#BD7848", "#C1654F", "#B14C33", "#913E2A", "#913E2A", "#803726", "#713122", "#6B2E20")
color_mapping <- colorRamp2(breakpoints, custom_colors)

png("heatmap_plot_5283_D80.png", width = 600, height = 1500)
Heatmap(plotdat[ordered_rows, ],  
        name = "Expression",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        show_column_names = FALSE,
        show_row_names = TRUE,
        bottom_annotation = col_ha,
        heatmap_legend_param = list(title = "Log\nNormalized\nExpression"),
        column_order = ordered_indices_in_plotdat,
        col = custom_colors,
        use_raster = FALSE,
        row_order = ordered_rows
)
dev.off()

