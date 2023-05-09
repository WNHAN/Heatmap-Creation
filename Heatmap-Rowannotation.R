library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
library(heatmaply)
library(ComplexHeatmap)
library(circlize)


#1. download OTU file from github human_16S.even13190.rel.sig.csv


##-- global variables
rare.depth = 13190

otu2 <- read.csv('https://raw.githubusercontent.com/cribioinfo/sci2017_analysis/master/data/human_16S.even13190.rel.sig.csv')

print(dim(otu2)) #63 42

sample2 <- read.csv('https://raw.githubusercontent.com/cribioinfo/sci2017_analysis/master/data/human_16S.sampleinfo.csv')

dim(sample2)

#2. file combination

library(magrittr)

otu2 <- otu2[-1] %>%
  t() %>%
  as.data.frame() %>%
  setNames((otu2[,1])) %>%
  tibble::rownames_to_column("Sample")

##-- file combination

sample3 <- sample2 %>%
  select(Sample, Response)

joined_data <- left_join(otu2,sample3, by = "Sample")

write.csv(joined_data, file="Heatmap Joined Data.csv",row.names=FALSE)

# import the saved csv
work_data <- read.csv(file="Heatmap Joined Data.csv",row.names = 1 )

#transpose work_data into 64/42(including response)
#final_data <- as.data.frame(t(work_data))

#final <- final_data[-64,]


#dim(final)
#str(final_data)
#head(final_data)

colnames(joined_data)[1] <- "" #delete the name of column 1:Sample
#colnames(joined_data)[65] <- "" #delet the name of column 65:Response


dim(joined_data)

#str(work_data)
#print(dim(work_data)) #42 64(63+response)  63/42

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
