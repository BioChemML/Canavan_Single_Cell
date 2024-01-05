library(Seurat)
library(readxl)
library(fishplot)

data <- read_excel('Fish_Plot_Final_4-1.xlsx')

timepoints <- c(0, 30, 80)
frac.table <- as.matrix(data[,-1])
parents <- rep(0, nrow(frac.table))

frac.table <- frac.table * 2000
fish <- createFishObject(frac.table, parents, timepoints = timepoints)
fish <- layoutClones(fish)

custom_colors <- c("#EF0000", "#FF6000", "#FFCF00", "#383bf5", "#50FFAF", "#BFFF40")
fish@col <- custom_colors

png("fish_plot.png", width = 1000, height = 600)

fishPlot(fish, shape = "spline", title.btm = "Cell Type Dynamics",
         cex.title = 0.5, vlines = timepoints, 
         vlab = c("Day 0", "Day 30", "Day 80"))

dev.off()
