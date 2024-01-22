library(Seurat)
library(reshape2)
library(ggplot2)

# load seurat objects for D80 WT and Can. These objects were created and processed in before.
WT <- readRDS("/location_where_you_input_file_is/D80_WT_1k.rds")
Can <- readRDS("/location_where_you_input_file_is/D80_Can_1k.rds")

# define our genes of interest for subsetting. 
genes <- c("Rap1gds1", "Atp1a2", "Plcb4", "Pdp1", "4930403P22Rik", "A730060N03Rik", 
           "Crybg1", "Adra2b", "Eps8l2", "Kcnj16", "9330182L06Rik", "Lama2", 
           "Shox2", "Pcp4", "Chchd10", "Igsf5", "Rbfox1")

# extract the counts matrix first, then subset the original count matrix to contain only our genes of interest. I do this step here to reduce the datasize for further analysis, so further data processing can be faster.
subset_counts_WT <- WT[["Spatial"]]@counts[genes, ]
subset_counts_Can <- Can[["Spatial"]]@counts[genes, ]

mean_exp_WT <- colMeans(subset_counts_WT)
mean_exp_Can <- colMeans(subset_counts_Can)

# create a dataframe based on the matrix. This step we are converting the data structure to a structure that can be used in ggplot2 for plotting heatmap.
mean_exp_df <- data.frame(Gene = genes, WT = mean_exp_WT, Can = mean_exp_Can)

# melt the dataframe
mean_exp_df_melted <- melt(mean_exp_df, id.vars = "Gene", variable.name = "Group", value.name = "MeanExpression")

mean_exp_df_melted$Gene <- factor(mean_exp_df_melted$Gene, levels = gene_order)

# assign the groups (WT/Can) to the cells for the annotations, so we can tell which cells are Can and which are WT when plotting the heatmap.
df_melted$Group <- ifelse(df_melted$Cell %in% colnames(sampled_WT), "WT", "Can")

# specify the desired order of showing the genes from top to bottom.
gene_order <- c("Rap1gds1", "Atp1a2", "Plcb4", "Pdp1", "4930403P22Rik", "A730060N03Rik", 
                "Crybg1", "Adra2b", "Eps8l2", "Kcnj16", "9330182L06Rik", "Lama2", 
                "Shox2", "Pcp4", "Chchd10", "Igsf5", "Rbfox1")

# set the Gene variable in the dataframe to this order
df_melted$Gene <- factor(df_melted$Gene, levels = gene_order)

# determine the midpoint for the color scale
mid_point <- max(mean_exp_df_melted$MeanExpression) / 2

# create the heatmap. The reason I choose these colors is because I used them when plotting the average heatmap. This is purely for consistency.
p <- ggplot(mean_exp_df_melted, aes(x = Group, y = Gene)) +
  geom_tile(aes(fill = MeanExpression), color = "white") +
  scale_fill_gradient2(low = "#5994CA", mid = "yellow", high = "red", midpoint = mid_point) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save the heatmap.
ggsave(output_path, p, width = 10, height = 6)
