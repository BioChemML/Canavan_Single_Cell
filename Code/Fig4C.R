ptm = Sys.time()
library(CellChat)
library(patchwork)

cellchat.WT <- readRDS("/location_where_you_input_file_is/D80_WT_1k.rds")
cellchat.Can <- readRDS("/location_where_you_input_file_is/D80_Can_1k.rds")

object.list <- list(WT = cellchat.WT, Can = cellchat.Can)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

cellchat

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
combined_plot <- gg1 + gg2
ggsave("total_nI_nS.png", plot = combined_plot, width = 10, height = 8, dpi = 300)

unique_clusters <- unique(unlist(cellchat@idents))
cluster_mapping <- read.csv("cluster_to_class_sorted_2.csv", stringsAsFactors = FALSE)

cluster_mapping$CellChatName <- gsub("^0*([0-9]+).*", "\\1", substr(cluster_mapping$Cluster, 1, 4))
cluster_mapping$CellChatName <- paste0("clst", cluster_mapping$CellChatName)

cluster_mapping$SimplifiedClass <- gsub("^([0-9]{2}).*", "\\1", cluster_mapping$Class)

mapping_subset <- cluster_mapping[cluster_mapping$CellChatName %in% unique_clusters, ]

object.list <- lapply(object.list, function(x) {
  object_clusters <- as.character(x@idents)
  object_mapping <- mapping_subset[mapping_subset$CellChatName %in% object_clusters, ]
  group.cellType <- factor(object_mapping$SimplifiedClass)
  mergeInteractions(x, group.cellType)
})

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))

png("merged_circle_Oligo_D80.png", width = 20*300, height = 16*300, res = 500)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, 
                   idents.use = "31", 
                   label.edge= T, edge.label.cex = 0.1, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

dev.off()

