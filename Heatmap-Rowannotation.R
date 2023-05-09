
library(ComplexHeatmap)
library(circlize)

# import the saved csv
work_data <- read.csv(file="Heatmap Joined Data.csv",row.names = 1 )

##-- convert to log scale 
otu_log = log10(round(work_data[1:63] * rare.depth) + 1)

##-- prep input for heatmap 
joined_data_plot = t(scale(t(as.matrix(otu_log)), scale=FALSE))

##-- add sample annotation 
joined_sample_anno = joined_data[,c('Response'),drop=FALSE]
joined_sample_anno = joined_sample_anno[order(match(row.names(joined_sample_anno),
                                      row.names(joined_data_plot))),,drop=F]
joined_plot_colors = c('NonResponder' = '#0000CC','Responder'='#CC0000')

joined_sample_anno_colors = list(
  Response = joined_plot_colors
)

plot_anno = rowAnnotation(df = joined_sample_anno,
                              col = joined_sample_anno_colors)

print(plot_anno)

##-- set up R plot display options in notebook
options(jupyter.plot_mimetypes = "image/png")
options(repr.plot.width = 10, repr.plot.height = 8)

##-- draw heatmap
row.title = paste0('Metastatic Melanoma (n=',nrow(joined_data_plot),') ')
col.title = paste0('OTUs (n=',ncol(joined_data_plot),') ')
myheatcol = colorRamp2(c(-0.8, 0, 0.8), c("blue", "white", "red"))

p1 = Heatmap(joined_data_plot,
             na_col = "#000000",
             col = myheatcol,
             rect_gp = gpar(col = '#000000'),
             show_heatmap_legend = TRUE,
             column_title = col.title,
             row_title = row.title,
             row_dend_width = unit(2, "cm"),
             column_dend_height = unit(2, "cm"),
             cluster_rows = TRUE,
             cluster_columns = TRUE,
             clustering_distance_rows = "euclidean",
             clustering_method_rows = "ward.D2",
             clustering_distance_columns = "euclidean",
             clustering_method_columns = "ward.D2",
             show_row_names = TRUE,
             column_names_side = 'bottom',
             column_names_gp = gpar(fontsize = 12),
             row_names_gp = gpar(fontsize = 8),
             left_annotation = plot_anno[1:42],
             heatmap_legend_param = 
               list(legend_direction = "horizontal",
                    legend_width = unit(5, "cm"), 
                    color_bar = 'continuous',
                    title = '')
)

draw(p1, annotation_legend_side = "left", heatmap_legend_side = "top")
