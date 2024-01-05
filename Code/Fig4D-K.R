D80_WT <- readRDS("/group/aawanggrp/siyulin/RDS_for_Analysis/D80_WT_1k.rds")
D80_Can <- readRDS("/group/aawanggrp/siyulin/RDS_for_Analysis/D80_Can_1k.rds")
D30_WT <- readRDS("/group/aawanggrp/siyulin/RDS_for_Analysis/D30_WT_1k.rds")
D30_Can <- readRDS("/group/aawanggrp/siyulin/RDS_for_Analysis/D30_Can_1k.rds")

met <- D80_WT@meta.data

df <- met %>% filter(matched_clusters == "clst5283")

p <- ggplot(met, aes(x = x, y = y)) + 
  geom_point(size=0.1, color = "#8F9495") + # size, color for the whole brain
  geom_point(data = df, aes(x = x, y = y), color = "red", size = 2) + # highlight the current cluster in red
  labs(title=paste('Highlighted cluster:', current_cluster)) +
  theme(axis.title = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=12)) +
  coord_fixed(ratio = 1)

ggsave("D80_WT.png", plot=p, width=15, height=15, dpi = 600)
