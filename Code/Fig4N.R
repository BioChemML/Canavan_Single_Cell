library(Seurat)
library(CellChat)
library(ggplot2) 
library(patchwork)
library(igraph)

WT_D80 <- readRDS("/group/aawanggrp/siyulin/Cell_Commute/D80_WT_1k.rds")
Can_D80 <- readRDS("/group/aawanggrp/siyulin/Cell_Commute/D80_Can_1k.rds")

WT_D80_interactions <- WT_D80@net$count["clst5268", ]
Can_D80_interactions <- Can_D80@net$count["clst5268", ]

# 2. Filter out clusters with zero interactions with clst5268
WT_D80_nonzero_clusters <- names(WT_D80_interactions[WT_D80_interactions > 0])
Can_D80_nonzero_clusters <- names(Can_D80_interactions[Can_D80_interactions > 0])

# 3. Identify clusters that are present in one group but not the other for interactions with clst5268
different_clusters_Can_D80 <- setdiff(Can_D80_nonzero_clusters, WT_D80_nonzero_clusters)
str(different_clusters_Can_D80)

different_clusters_WT_D80 <- setdiff(WT_D80_nonzero_clusters, Can_D80_nonzero_clusters)
str(different_clusters_WT_D80)

# 4. Combine the results
different_clusters <- union(different_clusters_WT_D80, different_clusters_Can_D80)
str(different_clusters)

# 5. Identify common clusters between WT_D80 and Can_D80
common_clusters <- intersect(WT_D80_nonzero_clusters, Can_D80_nonzero_clusters)

# 6. Among common clusters, identify those with different interaction counts with clst5268
different_counts_clusters_D80_outgoing <- common_clusters[WT_D80_interactions[common_clusters] != Can_D80_interactions[common_clusters]]

str(different_counts_clusters_D80_outgoing)

# 7. n_interaction for each category
unique_name_in_Can_D80_total <- sum(Can_D80_interactions[different_clusters_Can_D80])

unique_name_in_WT_D80_total <- sum(WT_D80_interactions[different_clusters_WT_D80])

sum_of_unique_total <- unique_name_in_Can_D80_total + unique_name_in_WT_D80_total

# Corrected 7. Same_name_but_different_number
same_name_diff_number_total_WT_D80 <- sum(WT_D80_interactions[different_counts_clusters_D80_outgoing])
same_name_diff_number_total_Can_D80 <- sum(Can_D80_interactions[different_counts_clusters_D80_outgoing])

list(
  Unique_name_in_Can_D80 = unique_name_in_Can_D80_total,
  Unique_name_in_WT_D80 = unique_name_in_WT_D80_total,
  Sum_of_unique = sum_of_unique_total,
  same_name_diff_number_total_WT_D80 = same_name_diff_number_total_WT_D80,
  same_name_diff_number_total_Can_D80 = same_name_diff_number_total_Can_D80
)

# Incoming
# 1. Extract Incoming Interaction Data for clst5268
WT_D80_incoming_interactions <- WT_D80@net$count[, "clst5268"]
Can_D80_incoming_interactions <- Can_D80@net$count[, "clst5268"]

# 2. Filter out clusters with zero incoming interactions with clst5268
WT_D80_nonzero_incoming_clusters <- names(WT_D80_incoming_interactions[WT_D80_incoming_interactions > 0])
Can_D80_nonzero_incoming_clusters <- names(Can_D80_incoming_interactions[Can_D80_incoming_interactions > 0])

# 3. Identify clusters that are present in one group but not the other for incoming interactions to clst5268
different_incoming_clusters_Can_D80 <- setdiff(Can_D80_nonzero_incoming_clusters, WT_D80_nonzero_incoming_clusters)
str(different_incoming_clusters_Can_D80)

different_incoming_clusters_WT_D80 <- setdiff(WT_D80_nonzero_incoming_clusters, Can_D80_nonzero_incoming_clusters)
str(different_incoming_clusters_WT_D80)

# 4. Combine the results
different_incoming_clusters <- union(different_incoming_clusters_WT_D80, different_incoming_clusters_Can_D80)
str(different_incoming_clusters)

# 5. Among common clusters, identify those with different incoming interaction counts to clst5268
common_incoming_clusters <- intersect(WT_D80_nonzero_incoming_clusters, Can_D80_nonzero_incoming_clusters)

different_incoming_counts_clusters_D80 <- common_incoming_clusters[WT_D80_incoming_interactions[common_incoming_clusters] != Can_D80_incoming_interactions[common_incoming_clusters]]

str(different_incoming_counts_clusters_D80)

# 8
unique_name_in_Can_D80_total <- sum(Can_D80_incoming_interactions[different_incoming_clusters_Can_D80])

unique_name_in_WT_D80_total <- sum(WT_D80_incoming_interactions[different_incoming_clusters_WT_D80])

sum_of_unique_total <- unique_name_in_Can_D80_total + unique_name_in_WT_D80_total

# Corrected 4. Same_name_but_different_number
same_name_diff_number_total_WT_D80 <- sum(WT_D80_incoming_interactions[different_incoming_counts_clusters_D80])
same_name_diff_number_total_Can_D80 <- sum(Can_D80_incoming_interactions[different_incoming_counts_clusters_D80])

list(
  Unique_name_in_Can_D80 = unique_name_in_Can_D80_total,
  Unique_name_in_WT_D80 = unique_name_in_WT_D80_total,
  Sum_of_unique = sum_of_unique_total,
  same_name_diff_number_total_WT_D80 = same_name_diff_number_total_WT_D80,
  same_name_diff_number_total_Can_D80 = same_name_diff_number_total_Can_D80
)

clusters <- c("clst1548", "clst181", "clst1837", "clst191", "clst26", "clst274", "clst330", "clst4808", "clst502", "clst5056", "clst550", "clst557", "clst578", "clst679", "clst816", "clst86", "clst5283")

# extract @netP
# Find the indices for clst1795 in senders
sender_names <- attr(Can_D80@netP$prob, "dimnames")[[1]]
receiver_names <- attr(Can_D80@netP$prob, "dimnames")[[2]]

sender_index <- which(sender_names %in% clusters)

# Find the indices for the clusters of interest in the receiver list
receiver_indices <- which(receiver_names %in% clusters)

# Subset the prob 3D array for interactions where clst1795 is the sender and the clusters in your list are the receivers
subset_prob <- Can_D80@netP$prob[sender_index, receiver_indices, , drop = FALSE]

# Check the structure of the subsetted prob to ensure it is a 3D array
str(subset_prob)

# Subset the prob array to only include the active pathways
final_prob <- subset_prob

# Subset the pathways
active_pathway_names <- Can_D80@netP$pathways

# Create the new netP structure
new_netP <- list(
  pathways = active_pathway_names,
  prob = final_prob
)

str(new_netP)


clusters <- c("clst5268", "clst5051", "clst5226")

# extract @idents
# Extract the indices of idents that match the clusters of interest
indices <- which(Can_D80@idents %in% clusters)

# Subset idents using these indices
subset_idents <- Can_D80@idents[indices]

# Drop unused levels
subset_idents <- factor(subset_idents)

# Check the structure of the subsetted idents
str(subset_idents)

# extract net
# net$prob
sender_names <- attr(Can_D80@net$prob, "dimnames")[[1]]
receiver_names <- attr(Can_D80@net$prob, "dimnames")[[2]]

sender_index <- which(sender_names %in% clusters)

# Find the indices for the clusters of interest in the receiver list
receiver_indices <- which(receiver_names %in% clusters)

net_prob_subset <- Can_D80@net$prob[sender_index, receiver_indices, , drop = FALSE]

str(net_prob_subset)

net_prob_final <- net_prob_subset

str(net_prob_final)

# net$pval
sender_names <- attr(Can_D80@net$pval, "dimnames")[[1]]
receiver_names <- attr(Can_D80@net$pval, "dimnames")[[2]]
sender_index <- which(sender_names %in% clusters)
receiver_indices <- which(receiver_names %in% clusters)
active_pathways <- dimnames(net_prob_final)[[3]]
net_pval_subset <- Can_D80@net$pval[sender_index, receiver_indices, active_pathways, drop = FALSE]
str(net_pval_subset)

# net$count
sender_names <- attr(Can_D80@net$count, "dimnames")[[1]]
receiver_names <- attr(Can_D80@net$count, "dimnames")[[2]]
sender_index <- which(sender_names %in% clusters)
receiver_indices <- which(receiver_names %in% clusters)
net_count_subset <- Can_D80@net$count[sender_index, receiver_indices, drop = FALSE]

# net#weight
sender_names <- attr(Can_D80@net$weight, "dimnames")[[1]]
receiver_names <- attr(Can_D80@net$weight, "dimnames")[[2]]
sender_index <- which(sender_names %in% clusters)
receiver_indices <- which(receiver_names %in% clusters)
net_weight_subset <- Can_D80@net$weight[sender_index, receiver_indices, drop = FALSE]

# put them together to make a new @net
new_net <- list(
  prob = net_prob_final,
  pval = net_pval_subset,
  count = net_count_subset,
  weight = net_weight_subset
)

str(new_net)

# extract LR
list3 <- dimnames(net_pval_subset)[[3]]
Can_D80_subset_LRsig <- Can_D80@LR$LRsig[Can_D80@LR$LRsig$interaction_name %in% list3, ]
new_LR <- list(LRsig = Can_D80_subset_LRsig)
str(new_LR)

# reconstruct the cellchat obj.
D80_WT_new <- new("CellChat")

# Assign the updated components to the new object
D80_WT_new@netP <- new_netP
D80_WT_new@idents <- subset_idents
D80_WT_new@net <- new_net
D80_WT_new@LR <- new_LR

D80_WT_new <- netAnalysis_computeCentrality(D80_WT_new, slot.name = "netP")
groupSize <- as.numeric(table(D80_WT_new@idents))

png("plot_Can.png", width = 40*300, height = 32*300, res = 500)
netVisual_circle(D80_WT_new@net$count, 
                 vertex.weight = groupSize,
                 idents.use = "clstclst5283",
                 title.name = "Number of interactions")
dev.off()

