library(Seurat)
library(CellChat)
library(ggplot2) 
library(patchwork)
library(igraph)

WT_D30 <- readRDS("/group/aawanggrp/siyulin/Cell_Commute/D30_WT_1k.rds")
Can_D30 <- readRDS("/group/aawanggrp/siyulin/Cell_Commute/D30_Can_1k.rds")

WT_D30_interactions <- WT_D30@net$count["clst5268", ]
Can_D30_interactions <- Can_D30@net$count["clst5268", ]

# 2. Filter out clusters with zero interactions with clst5268
WT_D30_nonzero_clusters <- names(WT_D30_interactions[WT_D30_interactions > 0])
Can_D30_nonzero_clusters <- names(Can_D30_interactions[Can_D30_interactions > 0])

# 3. Identify clusters that are present in one group but not the other for interactions with clst5268
different_clusters_Can_D30 <- setdiff(Can_D30_nonzero_clusters, WT_D30_nonzero_clusters)
str(different_clusters_Can_D30)

different_clusters_WT_D30 <- setdiff(WT_D30_nonzero_clusters, Can_D30_nonzero_clusters)
str(different_clusters_WT_D30)

# 4. Combine the results
different_clusters <- union(different_clusters_WT_D30, different_clusters_Can_D30)
str(different_clusters)

# 5. Identify common clusters between WT_D30 and Can_D30
common_clusters <- intersect(WT_D30_nonzero_clusters, Can_D30_nonzero_clusters)

# 6. Among common clusters, identify those with different interaction counts with clst5268
different_counts_clusters_D30_outgoing <- common_clusters[WT_D30_interactions[common_clusters] != Can_D30_interactions[common_clusters]]

str(different_counts_clusters_D30_outgoing)

# 7. n_interaction for each category
unique_name_in_Can_D30_total <- sum(Can_D30_interactions[different_clusters_Can_D30])

unique_name_in_WT_D30_total <- sum(WT_D30_interactions[different_clusters_WT_D30])

sum_of_unique_total <- unique_name_in_Can_D30_total + unique_name_in_WT_D30_total

# Corrected 7. Same_name_but_different_number
same_name_diff_number_total_WT_D30 <- sum(WT_D30_interactions[different_counts_clusters_D30_outgoing])
same_name_diff_number_total_Can_D30 <- sum(Can_D30_interactions[different_counts_clusters_D30_outgoing])


list(
  Unique_name_in_Can_D30 = unique_name_in_Can_D30_total,
  Unique_name_in_WT_D30 = unique_name_in_WT_D30_total,
  Sum_of_unique = sum_of_unique_total,
  same_name_diff_number_total_WT_D30 = same_name_diff_number_total_WT_D30,
  same_name_diff_number_total_Can_D30 = same_name_diff_number_total_Can_D30
)

# Incoming
# 1. Extract Incoming Interaction Data for clst5268
WT_D30_incoming_interactions <- WT_D30@net$count[, "clst5268"]
Can_D30_incoming_interactions <- Can_D30@net$count[, "clst5268"]

# 2. Filter out clusters with zero incoming interactions with clst5268
WT_D30_nonzero_incoming_clusters <- names(WT_D30_incoming_interactions[WT_D30_incoming_interactions > 0])
Can_D30_nonzero_incoming_clusters <- names(Can_D30_incoming_interactions[Can_D30_incoming_interactions > 0])

# 3. Identify clusters that are present in one group but not the other for incoming interactions to clst5268
different_incoming_clusters_Can_D30 <- setdiff(Can_D30_nonzero_incoming_clusters, WT_D30_nonzero_incoming_clusters)
str(different_incoming_clusters_Can_D30)

different_incoming_clusters_WT_D30 <- setdiff(WT_D30_nonzero_incoming_clusters, Can_D30_nonzero_incoming_clusters)
str(different_incoming_clusters_WT_D30)

# 4. Combine the results
different_incoming_clusters <- union(different_incoming_clusters_WT_D30, different_incoming_clusters_Can_D30)
str(different_incoming_clusters)

# 5. Among common clusters, identify those with different incoming interaction counts to clst5268
common_incoming_clusters <- intersect(WT_D30_nonzero_incoming_clusters, Can_D30_nonzero_incoming_clusters)

different_incoming_counts_clusters_D30 <- common_incoming_clusters[WT_D30_incoming_interactions[common_incoming_clusters] != Can_D30_incoming_interactions[common_incoming_clusters]]

str(different_incoming_counts_clusters_D30)

# 8
unique_name_in_Can_D30_total <- sum(Can_D30_incoming_interactions[different_incoming_clusters_Can_D30])

unique_name_in_WT_D30_total <- sum(WT_D30_incoming_interactions[different_incoming_clusters_WT_D30])

sum_of_unique_total <- unique_name_in_Can_D30_total + unique_name_in_WT_D30_total

# Corrected 4. Same_name_but_different_number
same_name_diff_number_total_WT_D30 <- sum(WT_D30_incoming_interactions[different_incoming_counts_clusters_D30])
same_name_diff_number_total_Can_D30 <- sum(Can_D30_incoming_interactions[different_incoming_counts_clusters_D30])



list(
  Unique_name_in_Can_D30 = unique_name_in_Can_D30_total,
  Unique_name_in_WT_D30 = unique_name_in_WT_D30_total,
  Sum_of_unique = sum_of_unique_total,
  same_name_diff_number_total_WT_D30 = same_name_diff_number_total_WT_D30,
  same_name_diff_number_total_Can_D30 = same_name_diff_number_total_Can_D30
)

clusters <- c("clst5267", "clst1008", "clst112", "clst1276", "clst1325", "clst1368", "clst1467", "clst2", "clst2485", "clst2785", "clst3799", "clst3946", "clst536")

# extract @netP
# Find the indices for clst1795 in senders
sender_names <- attr(Can_D30@netP$prob, "dimnames")[[1]]
receiver_names <- attr(Can_D30@netP$prob, "dimnames")[[2]]

sender_index <- which(sender_names %in% clusters)

# Find the indices for the clusters of interest in the receiver list
receiver_indices <- which(receiver_names %in% clusters)

# Subset the prob 3D array for interactions where clst1795 is the sender and the clusters in your list are the receivers
subset_prob <- Can_D30@netP$prob[sender_index, receiver_indices, , drop = FALSE]

# Check the structure of the subsetted prob to ensure it is a 3D array
str(subset_prob)

# Subset the prob array to only include the active pathways
final_prob <- subset_prob

# Subset the pathways
active_pathway_names <- Can_D30@netP$pathways

# Create the new netP structure
new_netP <- list(
  pathways = active_pathway_names,
  prob = final_prob
)

str(new_netP)

# extract @idents
# Extract the indices of idents that match the clusters of interest
indices <- which(Can_D30@idents %in% clusters)

# Subset idents using these indices
subset_idents <- Can_D30@idents[indices]

# Drop unused levels
subset_idents <- factor(subset_idents)

# Check the structure of the subsetted idents
str(subset_idents)

# extract net
# net$prob
sender_names <- attr(Can_D30@net$prob, "dimnames")[[1]]
receiver_names <- attr(Can_D30@net$prob, "dimnames")[[2]]

sender_index <- which(sender_names %in% clusters)

# Find the indices for the clusters of interest in the receiver list
receiver_indices <- which(receiver_names %in% clusters)

net_prob_subset <- Can_D30@net$prob[sender_index, receiver_indices, , drop = FALSE]

str(net_prob_subset)

net_prob_final <- net_prob_subset

str(net_prob_final)

# net$pval
sender_names <- attr(Can_D30@net$pval, "dimnames")[[1]]
receiver_names <- attr(Can_D30@net$pval, "dimnames")[[2]]
sender_index <- which(sender_names %in% clusters)
receiver_indices <- which(receiver_names %in% clusters)
active_pathways <- dimnames(net_prob_final)[[3]]
net_pval_subset <- Can_D30@net$pval[sender_index, receiver_indices, active_pathways, drop = FALSE]
str(net_pval_subset)

# net$count
sender_names <- attr(Can_D30@net$count, "dimnames")[[1]]
receiver_names <- attr(Can_D30@net$count, "dimnames")[[2]]
sender_index <- which(sender_names %in% clusters)
receiver_indices <- which(receiver_names %in% clusters)
net_count_subset <- Can_D30@net$count[sender_index, receiver_indices, drop = FALSE]

# net#weight
sender_names <- attr(Can_D30@net$weight, "dimnames")[[1]]
receiver_names <- attr(Can_D30@net$weight, "dimnames")[[2]]
sender_index <- which(sender_names %in% clusters)
receiver_indices <- which(receiver_names %in% clusters)
net_weight_subset <- Can_D30@net$weight[sender_index, receiver_indices, drop = FALSE]

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
Can_D30_subset_LRsig <- Can_D30@LR$LRsig[Can_D30@LR$LRsig$interaction_name %in% list3, ]
new_LR <- list(LRsig = Can_D30_subset_LRsig)
str(new_LR)

# reconstruct the cellchat obj.
D30_WT_new <- new("CellChat")

# Assign the updated components to the new object
D30_WT_new@netP <- new_netP
D30_WT_new@idents <- subset_idents
D30_WT_new@net <- new_net
D30_WT_new@LR <- new_LR

D30_WT_new <- netAnalysis_computeCentrality(D30_WT_new, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(D30_WT_new, width = 30, height = 30, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(D30_WT_new, width = 30, height = 30, pattern = "incoming")

png("plot_h1.png", width = 15*300, height = 17*300, res = 300)
print(ht1)
dev.off()

png("plot_h2.png", width = 15*300, height = 17*300, res = 300)
print(ht2)
dev.off()
