D80WT <- readRDS("/location_where_you_input_file_is/D80_WT_class.rds")
D80Can <- readRDS("/location_where_you_input_file_is/D80_Can_class.rds")

WT_subset <- subset(WT, subset = SimplifiedClass %in% c("23", "30", "31", "19"))
cell_ids_wt <- colnames(WT_subset)
WT_subset <- subset(WT, cells = cell_ids_wt)

Can_subset <- subset(Can, subset = SimplifiedClass %in% c("23", "30", "31", "19"))
cell_ids_can <- colnames(Can_subset)
Can_subset <- subset(Can, cells = cell_ids_wt)

WT_subset$orig.ident <- "D80_WT"
Can_subset$orig.ident <- "D80_Can"

D80 <- merge(WT_subset, y = Can_subset, add.cell.ids = c("D80WT", "D80Can"))

DefaultAssay(D80) <- 'Spatial'
D80 <- NormalizeData(D80, normalization.method = "LogNormalize", scale.factor = 10000)

D80 <- FindVariableFeatures(D80, selection.method = "vst", nfeatures = 8000)

all.genes <- rownames(D80)
D80 <- ScaleData(D80, features = all.genes)

D80 <- RunPCA(D80, dims = 1:30, features = VariableFeatures(object = D80))

D80 <- FindNeighbors(D80, dims = 1:30)
D80 <- RunUMAP(D80, dims = 1:30, res = 4)

SimplifiedClass <- D80@meta.data$SimplifiedClass
D80 <- SetIdent(D80, value = SimplifiedClass)

p2 <- DimPlot(D80, reduction = "umap", pt.size = 4, cols = colors, split.by = "orig.ident")
ggsave(filename = "D30_umap.png", plot = p2, width = 90, height = 50, limitsize = FALSE)
