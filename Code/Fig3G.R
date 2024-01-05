library(rbioapi)

gene_new <- c("Ndufa10", "Dlst", "Ndufs1", "Hspd1", "Serpinb7", "Il22ra2", "Casp8", "Gpx3", "Il10", "Il1r1", "Il1r2", "Trpa1", "Creb1", "Irs1", "Cd28", "Icos", "Ptprc", "Batf3", "Trpa1", "Arhgap30", "Ifngr1", "IL18r1")

proteins_mapped <- rba_string_map_ids(ids = gene_new, 
                                      species = 10090)
protein_use <- proteins_mapped$stringId

graph_1 <- rba_string_network_image(ids = protein_use,
                                    image_format = "highres_image",
                                    species = 10090,
                                    required_score = 500,
                                    network_flavor = "confidence",
                                    hide_node_labels = FALSE,
                                    hide_disconnected_nodes = TRUE,
                                    hide_structure_pics = TRUE)