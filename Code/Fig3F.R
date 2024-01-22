library(Seurat)
library(reshape2)
library(ggplot2)

WT <- readRDS("location_where_you_input_file_is/D30_WT_1k.rds")
Can <- readRDS("location_where_you_input_file_is/D30_Can_1k.rds")

WT$orig.ident <- "D30_WT"
Can$orig.ident <- "D30_Can"

D30 <- merge(WT, y = Can)

D30 <- NormalizeData(D30, normalization.method = "LogNormalize", scale.factor = 10000)

D30 <- FindVariableFeatures(D30, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(D30)
D30 <- ScaleData(D30, features = all.genes)

D30 <- RunPCA(D30, features = VariableFeatures(object = D30))
print(D30[["pca"]], dims = 1:5, nfeatures = 5)

de.results_D30_2 <- FindMarkers(D30,
                                ident.1 = "D30_Can", 
                                ident.2 = "D30_WT",
                                group.by = "orig.ident",
                                test.use = "MAST",
                                min.pct = 0.001,
                                logfc.threshold = 0.25)

write.xlsx(de.results_D30_2, file = "Fig3F_all.xlsx", rowNames = TRUE)

genes_of_interest <- c("Rgs20", "Wdr12", "Olfr1413", "Cflar", "Rgs16", "Gpx3", "Ggt1", "Erbb4", "Ifngr1", "Pik3r5", "Ndufs1", "Hspd1", "Rgs5", "Cd28", "Irs1", "Il10", "Tmeff2", "Serpinb7", "Olfr218", "Matk", "Gng7", "Ndufa10", "Pak1ip1", "Olfr417", "Aim2", "Creb1", "Dlst", "Il22ra2", "Arhgap30", "Casp8", "Batf3", "Ptprc", "Il1r1", "Imp4", "Trpa1", "Il18r1", "Icos", "Prkag3")

scaled_data <- GetAssayData(D30, slot = "scale.data")[genes_of_interest, ]
D30$Condition <- ifelse(D30$orig.ident == "D30_Can", "D30_Can", "D30_WT")
cells_can <- colnames(scaled_data)[D30$Condition == "D30_Can"]
cells_wt <- colnames(scaled_data)[D30$Condition == "D30_WT"]
set.seed(200) 
sample_cells_can <- sample(cells_can, size = round(length(cells_can) * 1))
sample_cells_wt <- sample(cells_wt, size = round(length(cells_wt) * 1))
sample_cells <- c(sample_cells_can, sample_cells_wt)
conditions <- c(rep("D30_Can", length(sample_cells_can)), rep("D30_WT", length(sample_cells_wt)))
names(conditions) <- sample_cells
sample_cells_sorted <- sample_cells[order(conditions)]
plotdat <- as.matrix(scaled_data[, sample_cells_sorted])


outlier_cells <- colnames(plotdat)[apply(plotdat, 2, function(cell_data) any(cell_data > 100))]
plotdat <- plotdat[, !(colnames(plotdat) %in% outlier_cells)]
conditions <- conditions[!(names(conditions) %in% outlier_cells)]

col_ha <- HeatmapAnnotation(df = data.frame(Condition = conditions[sample_cells_sorted]))
can_indices_in_plotdat <- which(D30$orig.ident[colnames(plotdat)] == "D30_Can")
wt_indices_in_plotdat <- which(D30$orig.ident[colnames(plotdat)] == "D30_WT")
ordered_indices_in_plotdat <- c(can_indices_in_plotdat, wt_indices_in_plotdat)

gene_order <- c("Rgs20", "Wdr12", "Olfr1413", "Cflar", "Rgs16", "Gpx3", "Ggt1", "Erbb4", "Ifngr1", "Pik3r5", "Ndufs1", "Hspd1", "Rgs5", "Cd28", "Irs1", "Il10", "Tmeff2", "Serpinb7", "Olfr218", "Matk", "Gng7", "Ndufa10", "Pak1ip1", "Olfr417", "Aim2", "Creb1", "Dlst", "Il22ra2", "Arhgap30", "Casp8", "Batf3", "Ptprc", "Il1r1", "Imp4", "Trpa1", "Il18r1", "Icos", "Prkag3")

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

png("heatmap_plot_C30_D80_4.png", width = 600, height = 800)
Heatmap(plotdat[ordered_rows, ],  
        name = "Expression",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        show_column_names = FALSE,
        show_row_names = TRUE,
        bottom_annotation = col_ha,
        heatmap_legend_param = list(title = "Log\nNormalized\nExpression"),
        column_order = ordered_indices_in_plotdat,
        col = color_mapping,
        use_raster = FALSE,
        row_order = ordered_rows
)
dev.off()
