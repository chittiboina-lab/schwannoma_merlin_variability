#### Plotting the study summary for VS Schwannoma paper ####
#### Install libraries ####
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("catmaply")
install.packages("circlize")
install.packages("ComplexHeatmap")

#### Load packages and set seed ####
library(ggplot2)
library(tidyverse)
library(catmaply)
library(pheatmap)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
set.seed(42)

#### Load data ####
df <- read.csv("data/metadata_vs_study_updated.csv",header = TRUE, na.strings = "-")
rownames(df) <- df$Sample
df_small <- df[c("Sample", "Patient", "Side", "Diagnosis", "Sex", "Prep",
                 "scRNAseq", "snRNAseq", "Multiome",
                 "RNAScope", "mIHC", "spatialTx", "Methylation",
                 "Growth_rate", "Tumor_volume", 
                 "Genotype" )]
str(df_small)
df_small$Patient <- as.character(df_small$Patient)
mydata <- as_tibble(df_small)
mydata
data <- mydata %>% select(scRNAseq:Methylation)
data
data <- t(data)

#### Create top annotation ####
# # Dynamically create colors for patients
# nb.cols <- length(unique(df_small$Patient))
# my.colors <- colorRampPalette(brewer.pal(min(nb.cols, 9),'Set1'))(nb.cols)
# names(my.colors) <- unique(df_small$Patient)

my.colors <- c( "1" = '#e74c3c', "2" = '#d35400', "3" =  '#f39c12',
"4" = '#cd6155', "5" =  '#a93226', "6" =  '#ba4a00', "7" = '#f1c40f', "8" = '#d68910', "9" = '#b9770e', 
"10" = '#3498db', "11" = '#2980b9', "12" =  '#16a085',
"13" = '#8e44ad', "14" = '#34495e', "15" = '#5dade2')

# Update col_top list using my.colors for Patient
col_top = list(
  Patient = my.colors,
  Side = c("Right" = "#FF0000", "Left" = "#FFFF00"),
  Diagnosis = c('NF2'= '#e74c3c', 'Sporadic'= '#3498db'),
  Sex = c("M" = "#66FFFF", "F" = "#FF66B2"),
  Prep = c("fresh + archival" = "#FFFF66", "fresh" = "#66FFFF", "archival" = "#FF6666")
)

# HeatmapAnnotation with the updated col_top
ha_top <- HeatmapAnnotation(
  Patient = df_small$Patient,
  Side = df_small$Side,
  Diagnosis = df_small$Diagnosis,
  Sex = df_small$Sex,
  Prep = df_small$Prep,
  border = TRUE, 
  col = col_top,
  gap = unit(1, "mm"),
  gp = gpar(col = "white", lwd = 2),
  show_legend = c(Patient = FALSE, Side = TRUE, Diagnosis = TRUE, Sex = TRUE, Prep = TRUE)  # Suppress legend for Patient
)

#### Create bottom annotation ####
col_btm = list( Growth.rate = circlize::colorRamp2(c(3, 20, 57), 
                                                 c("deepskyblue", "#FFFF99", "#FF9933")),
                Tumor.volume = circlize::colorRamp2(c(2400, 10000, 25090), 
                                              c("deepskyblue", "#FFFF99", "#FF9933")),
                Genotype = c("not_found" = "black","nonsense/truncating" = "red",
                "splicing" = "blue", "NA" = "white")
)


ha_btm <- HeatmapAnnotation(  Growth.rate = df_small$Growth.rate,
                              Tumor.volume = df_small$Tumor.volume,
                              Genotype = df_small$Genotype,
                              col = col_btm,
                              na_col = "white", border = TRUE,
                              gap = unit(1, "mm"),
                              gp = gpar(col = "white", lwd = 2)
)


# Create a named vector for color mapping
color_mapping <- c("Y" = "#009900", "NA" = "white")

#### Plot the integrated heatmap and save ####
p1 <- Heatmap(
  data,
  col = color_mapping,  # Use the named vector for color mapping
  na_col = "white",
  top_annotation = ha_top,
  height = unit(2.5, "cm"),
  bottom_annotation = ha_btm,
  show_heatmap_legend = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",
  show_row_names = TRUE,
  column_labels = df$Sample,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 2)
)


pdf("plots/test_13.pdf", height = 5, width = 11, fonts = "Helvetica")
p1
dev.off()

#### Session info ####
sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.4
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.18.0 RColorBrewer_1.1-3    circlize_0.4.15      
# [4] pheatmap_1.0.12       catmaply_0.9.3        lubridate_1.9.3      
# [7] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4          
# [10] purrr_1.0.2           readr_2.1.4           tidyr_1.3.0          
# [13] tibble_3.2.1          tidyverse_2.0.0       ggplot2_3.4.4        
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.4          generics_0.1.3      shape_1.4.6         stringi_1.8.2      
# [5] hms_1.1.3           digest_0.6.33       magrittr_2.0.3      timechange_0.2.0   
# [9] iterators_1.0.14    foreach_1.5.2       doParallel_1.0.17   GlobalOptions_0.1.2
# [13] BiocManager_1.30.22 fansi_1.0.6         scales_1.3.0        codetools_0.2-19   
# [17] cli_3.6.1           crayon_1.5.2        rlang_1.1.2         munsell_0.5.0      
# [21] withr_2.5.2         tools_4.3.1         parallel_4.3.1      tzdb_0.4.0         
# [25] colorspace_2.1-0    BiocGenerics_0.48.1 GetoptLong_1.0.5    vctrs_0.6.5        
# [29] R6_2.5.1            png_0.1-8           stats4_4.3.1        matrixStats_1.1.0  
# [33] lifecycle_1.0.4     S4Vectors_0.40.2    IRanges_2.36.0      clue_0.3-65        
# [37] cluster_2.1.4       pkgconfig_2.0.3     pillar_1.9.0        gtable_0.3.4       
# [41] glue_1.6.2          tidyselect_1.2.0    rstudioapi_0.15.0   rjson_0.2.21       
# [45] compiler_4.3.1  
# 
