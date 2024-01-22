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

# randomly select 30 cells from each count matrixs
sampled_WT <- subset_counts_WT[, sample(1:ncol(subset_counts_WT), 30)]
sampled_Can <- subset_counts_Can[, sample(1:ncol(subset_counts_Can), 30)]

# combine the selected count matrixs(up to this point, the dataset we are going to use will have 60 cells and 17 genes. 30 of the cells from WT and 30 from Can)
combined_sampled_data <- as.matrix(cbind(sampled_WT, sampled_Can))

# create a dataframe based on the matrix. This step we are converting the data structure to a structure that can be used in ggplot2 for plotting heatmap.
df <- as.data.frame(as.matrix(combined_sampled_data))

# this step is to give our dataframe's rowname a name. After this conversion, the rownames will be considered as a "column", so it won't be lost during the melting process.
df$Gene <- rownames(df)

# melt the dataframe
df_melted <- melt(df, id.vars = "Gene", variable.name = "Cell", value.name = "Expression")

# check the top rows of the melted dataframe. I did this to make sure the datastructure is what I want.
head(df_melted)

# assign the groups (WT/Can) to the cells for the annotations, so we can tell which cells are Can and which are WT when plotting the heatmap.
df_melted$Group <- ifelse(df_melted$Cell %in% colnames(sampled_WT), "WT", "Can")

# specify the desired order of showing the genes from top to bottom.
gene_order <- c("Rap1gds1", "Atp1a2", "Plcb4", "Pdp1", "4930403P22Rik", "A730060N03Rik", 
                "Crybg1", "Adra2b", "Eps8l2", "Kcnj16", "9330182L06Rik", "Lama2", 
                "Shox2", "Pcp4", "Chchd10", "Igsf5", "Rbfox1")

# set the Gene variable in the dataframe to this order
df_melted$Gene <- factor(df_melted$Gene, levels = gene_order)

# determine the midpoint for the color scale
mid_point <- max(df_melted$Expression) / 2

# define the path and file name where you want to save the plot.
output_path <- "location_where_you_output_file_is/single_cell_heatmap.png"

# create the heatmap. The reason I choose these colors is because I used them when plotting the average heatmap. This is purely for consistency.
p <- ggplot(df_melted, aes(x = Cell, y = Gene)) +
  geom_tile(aes(fill = Expression), color = "white") +
  scale_fill_gradient2(low = "#5994CA", mid = "yellow", high = "red", midpoint = mid_point) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(~Group, scales = "free", space = "free")

# save the heatmap.
ggsave(output_path, p, width = 10, height = 6)
