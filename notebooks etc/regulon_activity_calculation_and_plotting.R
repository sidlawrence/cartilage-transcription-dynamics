#libraries
library(dorothea)
library(ComplexHeatmap)
library(dplyr)
library(viper)
library(tibble)
library(tidyr)
library(circlize)


cart<-readRDS("path/to/object.rds")
cart<-SetIdent(cart, value = "annotation")
regulon <- get(data("dorothea_hs", package = "dorothea"))
cart <- run_viper(cart, regulon,options = list(method = "scale", minsize = 4, eset.filter = FALSE, cores = 1,   verbose = T))
DefaultAssay(object = cart) <- "dorothea"
cart <- ScaleData(cart)

viper_scores_df <- GetAssayData(cart, assay = "dorothea", slot = 'scale.data') %>% data.frame(check.names = F) %>% t()
CellsClusters <- data.frame(cell = names(Idents(cart)), cell_type = as.character(Idents(cart)),  check.names = F)
viper_scores_clusters <- viper_scores_df  %>% data.frame() %>%  rownames_to_column("cell") %>% gather(tf, activity, -cell) %>% inner_join(CellsClusters)

#summarize scores by cell annotation
summarized_viper_scores <- viper_scores_clusters %>% group_by(tf, cell_type) %>% summarise(avg = mean(activity), std = sd(activity))

#format for ease of plotting
highly_variable_tfs <- summarized_viper_scores %>% group_by(tf) %>%  mutate(var = var(avg))  %>%  ungroup() %>%  distinct(tf)
summarized_viper_scores_df <- summarized_viper_scores %>%  semi_join(highly_variable_tfs, by = "tf") %>%  dplyr::select(-std) %>%  spread(tf, avg) %>%  data.frame(row.names = 1, check.names = FALSE) 
summarized_viper_scores_df_t<- t(summarized_viper_scores_df)

#select genes for plotting #eg KLF1, SOX9, COL2A1, COL9A1
genes<-c("KLF1", "SOX9", "COL2A1", "COL9A1")
plotting_viper_scores_df<-summarized_viper_scores_df_t[(which(rownames(summarized_viper_scores_df_t) %in% genes)),]

#reorder as per list of genes
#order gene names for plotting
plotting_viper_scores_df<- plotting_viper_scores_df[match(genes, rownames(plotting_viper_scores_df)),]

#plot
Heatmap(plotting_viper_scores_df,col = col_fun, cluster_columns = F,cluster_rows = F, row_names_gp = gpar(fontsize = 6), rect_gp = gpar(col = "white", lwd = 0.2), column_title = "Regulon Activity", column_gap = unit(3, "mm"),show_row_names = T,cluster_column_slices = T, border=T, column_title_gp = gpar(fontsize = 20, fontface = "bold"))


