library(readxl)
library(ggplot2)
library(pheatmap)

data <- read_excel("C:/All data/Research/Canavan/StereoMap/Figure2_No_Nucleus.xlsx")

data_long <- reshape2::melt(data)

# Determine the midpoint value, for example:
midpoint_ori <- median(data_long$value, na.rm = TRUE)

# Create heatmap
Heatmap_2 <- ggplot(data_long, aes(x = variable, y = Gene, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "royalblue4", high = "red", mid = "white", midpoint = midpoint_ori+20, limit = c(min(data_long$value), max(data_long$value))) + 
  theme_minimal() +
  xlab("Column Name") +
  ylab("Row Name")

ggsave(filename = "C:/All data/Research/Canavan/StereoMap/plots/Figure2_Heatmap_No_Nucleus.png", plot = Heatmap_2, width = 4, height = 16, limitsize = FALSE)

