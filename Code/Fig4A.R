library(Seurat)
library(dplyr)
library(ggplot2)

can <- readRDS("/group/aawanggrp/siyulin/RDS_for_Analysis/D80_Can_1k.rds") # or load WT
wt <- readRDS("/group/aawanggrp/siyulin/RDS_for_Analysis/D80_WT_1k.rds")

# extract the meta.data. It has the cluster identities stored here.
met <- can@meta.data

# extract the cluster column.
all_clusters <- unique(met$matched_clusters)

# define the clusters we are going to show. I selected based on the region shown in the excel. They are
# all good clusters that are presented in one region only.
highlight_clusters <- c("clst1795", "clst3191", "clst5267", "clst232", "clst233", "clst237", "clst184", 
                        "clst188", "clst190", "clst200", "clst177", "clst182", "clst174", "clst176", 
                        "clst82", "clst99", "clst69", "clst70", "clst72", "clst75", "clst140", "clst65", 
                        "clst66", "clst17", "clst18", "clst23", "clst27", "clst30", "clst49", "clst59", 
                        "clst328", "clst329", "clst433", "clst420", "clst421", "clst940", "clst1062", 
                        "clst1072", "clst1075", "clst1088", "clst1089", "clst1090", "clst1092", "clst1093", 
                        "clst418", "clst1403", "clst1405", "clst1407", "clst2878", "clst2880", "clst2881", 
                        "clst2882", "clst2884", "clst2888", "clst2890", "clst2891", "clst2893", "clst2894", 
                        "clst2904", "clst2907", "clst2909", "clst1056", "clst1069", "clst1191", "clst1912", 
                        "clst1932", "clst2465", "clst2467", "clst1792", "clst1793", "clst1794", "clst3155", 
                        "clst5230", "clst5231", "clst5232", "clst5233", "clst314", "clst413", "clst3374", 
                        "clst5021", "clst5022", "clst5023", "clst5032", "clst5033", "clst5039", "clst5040", 
                        "clst5041", "clst5043", "clst5044", "clst5045", "clst5049", "clst5050", "clst5051", 
                        "clst5058", "clst5061", "clst5065", "clst5066", "clst5067", "clst5068", "clst5075", 
                        "clst5076", "clst5088", "clst5090", "clst5100", "clst5171", "clst5184", "clst5190", 
                        "clst5192", "clst108", "clst139", "clst140", "clst141", "clst167", "clst168", "clst169",
                        "clst192", "clst205", "clst214", "clst228", "clst231", "clst233", "clst288", "clst393",
                        "clst400", "clst403", "clst430", "clst432", "clst434", "clst599", "clst134", "clst401",
                        "clst406", "clst969", "clst80", "clst120", "clst123", "clst124", "clst143", "clst187",
                        "clst5179", "clst4098")

# This step is just for in case. It is just for plotting some clusters that we want to further highlight. But I didn't
# use it when making the spatial plot.
big_clusters <- "clst1795"

# Extract the cell spatial data only about the clusters we defined earlier.
df <- met %>% filter(matched_clusters %in% highlight_clusters)

# make the spatial plot.
ggplot(met, aes(x = x, y = y)) + 
  geom_point(size=0.1, color = "#E2D67D") + # size, color for the whole brain
  geom_point(data = subset(df, !(matched_clusters %in% big_clusters)), aes(x = x, y = y, colour = matched_clusters), size = 0.1) +#for the clusters in the matched_clusters but not big_clusters, define their sizes.
  geom_point(data = subset(df, matched_clusters %in% big_clusters), aes(x = x, y = y, colour = matched_clusters), size = 0.1) +#for the clusters in the big_clusters but not matched_clusters, define their sizes.
  scale_color_manual(values = c("clst232" = "#9F0000", "clst233" = "#9F0000", "clst237" = "#9F0000", #define the colors for each cluster. I used the same color for the clusters in one region.
                                "clst205" = "#9F0000", "clst214" = "#9F0000", "clst228" = "#9F0000",
                                "clst231" = "#9F0000", "clst233" = "#9F0000",
                                "clst184" = "#9F0000", "clst190" = "#9F0000", "clst187" = "#9F0000",
                                "clst200" = "#9F0000", "clst177" = "#9F0000", 
                                "clst182" = "#E57373", "clst188" = "#E57373", "clst393" = "#E57373",
                                "clst174" = "#E57373", "clst176" = "#E57373", "clst167" = "#E57373",
                                "clst168" = "#E57373", "clst169" = "#E57373", "clst192" = "#E57373",
                                "clst400" = "#E57373", "clst403" = "#E57373", "clst599" = "#E57373",
                                "clst401" = "#E57373", "clst406" = "#E57373", "clst80" = "#E57373",
                                "clst4098" = "#E57373",
                                "clst82" = "#E57373", "clst99" = "#E57373", 
                                "clst69" = "#F7DA00", "clst70" = "#F7DA00", "clst72" = "#F7DA00", 
                                "clst75" = "#F7DA00", "clst140" = "#F7DA00", "clst65" = "#F7DA00", 
                                "clst66" = "#F7DA00", "clst17" = "#F7DA00", "clst18" = "#F7DA00", 
                                "clst23" = "#F7DA00", "clst27" = "#F7DA00", "clst30" = "#F7DA00", 
                                "clst49" = "#F7DA00", "clst59" = "#F7DA00", "clst108" = "#F7DA00", 
                                "clst139" = "#F7DA00", "clst140" = "#F7DA00", "clst141" = "#F7DA00",
                                "clst134" = "#F7DA00", "clst120" = "#F7DA00", "clst123" = "#F7DA00",
                                "clst124" = "#F7DA00", "clst143" = "#F7DA00",
                                "clst328" = "#67CEFF", "clst329" = "#67CEFF", "clst434" = "#67CEFF",
                                "clst433" = "#9AEDFF", "clst420" = "#B5FFFF", "clst421" = "#B5FFFF", 
                                "clst940" = "#7B68EE", "clst1062" = "#7B68EE", "clst1072" = "#7B68EE", 
                                "clst1075" = "#7B68EE", "clst1088" = "#7B68EE", "clst1089" = "#7B68EE", 
                                "clst1090" = "#7B68EE", "clst1092" = "#7B68EE", "clst1093" = "#7B68EE",
                                "clst288" = "#7B68EE", "clst969" = "#7B68EE",
                                "clst418" = "#FFD300", "clst1403" = "#FFD300", "clst1405" = "#FFD300", 
                                "clst1407" = "#FFD300", "clst430" = "#FFD300", "clst432" = "#FFD300",
                                "clst2878" = "#68CFB0", "clst2880" = "#68CFB0", "clst2881" = "#68CFB0", 
                                "clst2882" = "#68CFB0", "clst2884" = "#68CFB0", "clst2888" = "#68CFB0", 
                                "clst2890" = "#68CFB0", "clst2891" = "#68CFB0", "clst2893" = "#68CFB0", 
                                "clst2894" = "#68CFB0", "clst2904" = "#68CFB0", "clst2907" = "#68CFB0", 
                                "clst2909" = "#68CFB0",
                                "clst1056" = "#BDFEE0", "clst1069" = "#BDFEE0", "clst1191" = "#BDFEE0", 
                                "clst1912" = "#BDFEE0", "clst1932" = "#BDFEE0", "clst2465" = "#BDFEE0", 
                                "clst2467" = "#BDFEE0", 
                                "clst1792" = "#9400D3", "clst1793" = "#9400D3", "clst1794" = "#9400D3", 
                                "clst1797" = "#9400D3", "clst3155" = "#9400D3", "clst5230" = "#9400D3", 
                                "clst5231" = "#9400D3", "clst5232" = "#9400D3", "clst5233" = "#9400D3",
                                "clst5179" = "#9400D3", "clst1795" = "#9400D3",
                                "clst314" = "grey", "clst413" = "grey", 
                                "clst3374" = "#FF00FF", 
                                "clst5021" = "#FF00FF", "clst5022" = "#FF00FF", "clst5023" = "#FF00FF", 
                                "clst5032" = "#FF00FF", "clst5033" = "#FF00FF", "clst5039" = "#FF00FF", 
                                "clst5040" = "#FF00FF", "clst5041" = "#FF00FF", "clst5043" = "#FF00FF", 
                                "clst5044" = "#FF00FF", "clst5045" = "#FF00FF", "clst5049" = "#FF00FF", 
                                "clst5050" = "#FF00FF", "clst5051" = "#FF00FF", "clst5058" = "#FF00FF", 
                                "clst5061" = "#FF00FF", "clst5065" = "#FF00FF", "clst5066" = "#FF00FF", 
                                "clst5067" = "#FF00FF", "clst5068" = "#FF00FF", "clst5075" = "#FF00FF", 
                                "clst5076" = "#FF00FF", "clst5088" = "#FF00FF", "clst5090" = "#FF00FF", 
                                "clst5100" = "#FF00FF", "clst5171" = "#FF00FF", "clst3191" = "#FF00FF",
                                "clst5184" = "#6B4522", "clst5190" = "#6B4522", "clst5192" = "#6B4522" 
  ) )+ 
  labs(title='highlighted clusters') +
  theme(
    panel.background = element_rect(fill = "black"),# make the background black, because black background looks better and is consistent with our other plots.
    plot.background = element_rect(fill = "black"),# same as above
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, color = "white"), # text size here
    axis.title = element_blank(),
    plot.title = element_text(color = "white")
  ) +
  
  coord_fixed(ratio = 1)  # ensure a 1:1 aspect ratio. If you don't define it here, the spatial plot might be squeezed.

ggsave("/group/aawanggrp/siyulin/Plots/Spatial_Plots/Can_final.png", width = 15, height = 15, dpi = 900)