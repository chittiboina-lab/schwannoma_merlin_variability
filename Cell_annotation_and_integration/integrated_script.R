## SoupX test 
#setwd('~/OneDrive - National Institutes of Health/Projects/VS_integrated_0123/integrated_NF2_vs_Sporadic/')

{library(Seurat)
  library(SoupX)
  library(DropletUtils)
  library(ggplot2)
  library(DoubletFinder)
  library(Matrix)
  library(knitr)
  library(scSorter)
}
dir()
names <- c()
#needs operator with Gene Expression
paths <- c('snVS01.h5','snVS02.h5','snVS03.h5','snVS04.h5','snVS05.h5',
           'snVS06.h5','snVS07.h5','snVS08.h5','snVS09.h5')

for (path in paths){
  name <- substr(path,1,6)
  filt.matrix <- Read10X_h5(path,use.names = T)
    filt.matrix <- filt.matrix$`Gene Expression`
  names <- c(names, name)
  
  raw.matrix  <- Read10X_h5(stringr::str_interp("${name}_raw.h5"),use.names = T)
   raw.matrix <- raw.matrix$`Gene Expression`
  #str(raw.matrix)
  #str(filt.matrix)
  
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  #soup.channel
  
  b <- CreateSeuratObject(counts = filt.matrix, assay = "RNA", project = stringr::str_interp("${name}"))
  b <- NormalizeData(b)
  b <- FindVariableFeatures(b)
  b <- ScaleData(b)
  b <- RunPCA(b)
  b <- RunUMAP(b, dims = 1:30)
  b <- FindNeighbors(b, dims = 1:30)
  b <- FindClusters(b, resolution = 0.5)
  meta <- b@meta.data
  umap <- b@reductions$umap@cell.embeddings
  rm(b)
  
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  soup.channel = autoEstCont(soup.channel)
  #ggsave(stringr::str_interp('${name}/${name}_soup_channel.png'), plot = autoEstCont(soup.channel))
  
  adjusted_counts <- adjustCounts(soup.channel)
  obj <- CreateSeuratObject(counts = adjusted_counts, assay = "RNA", project = stringr::str_interp("${name}"))
  obj$type ="NF2"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('${name}/${name}_QC_stats.png'), plot = QC_stats)
  #QC_stats
  ## Remove empty droplets and low quality cells
  obj <- subset(obj, nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 25)
  post_cutoff_QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('${name}/${name}_post_cutoff_QC_stats.png'), plot = post_cutoff_QC_stats)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  ## Scale
  obj <- ScaleData(obj, features = rownames(data))
  ## Run PCA
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = 1:30)
  ## Find and Remove Doublets
  # pK Identification (no ground-truth)
  sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  #ggsave(stringr::str_interp('${name}/${name}_pK_identify.png'), plot = bcmvn_obj)
  annotations <- Idents(obj)
  ## ex: annotations <- obj@meta.data$ClusteringResults
  
  homotypic.prop <- modelHomotypic(annotations)         
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  
  nExp_poi <- round(0.075*nrow(obj@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder
  obj <- doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) ## Unsure of most appropriate parameters, these are from the example on the doublet finder github
  # this column in meta.data have a different name for every sample and/or parameter setting, so I store the name of the column to refer to it later without hard coding
  
  doublet_col <- rev(names(obj@meta.data))[1] 
  doublets <- DimPlot(obj, reduction = 'umap', group.by =doublet_col)
  ggsave(stringr::str_interp('${name}/${name}_doublets.png'), plot = doublets)
  
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5)
  ## Visualize before analysis
  clusters <- DimPlot(obj, reduction = 'umap')
  
  ggsave(stringr::str_interp('${name}/${name}_clusters.png'), plot = clusters)
  
  markerlist <- FindAllMarkers(obj)
  write.csv(markerlist, file = stringr::str_interp("${name}/top_markerlist.csv"))
  
  ## Using scSorter
  topgenes <- head(VariableFeatures(obj), 2000)
  obj_exp = GetAssayData(obj)
  topgene_filter = rowSums(as.matrix(obj_exp)[topgenes, ]!=0) > ncol(obj_exp)*.1
  topgenes = topgenes[topgene_filter]
  ## load annotation file
  rm(anno)
  anno <- read.csv("sorted_sc_markers_DA.csv")
  anno
  picked_genes = unique(c(anno$Marker, topgenes))
  write.csv(picked_genes,file = stringr::str_interp("${name}/${name}_picked_genes.csv"))
  obj_exp = obj_exp[rownames(obj_exp) %in% picked_genes, ]
  ## run scSorter
  rts <- scSorter(obj_exp, anno)
  rts
  obj$pred.type <- rts$Pred_Type
  table(obj$pred.type)
  
  cell_types <- DimPlot(obj,  reduction = 'umap' , group.by = "pred.type")
  ggsave(stringr::str_interp('${name}/${name}_celltypes.png'), plot = cell_types)
  ## Tabulate Cell Types
  write.csv(table(Idents(obj)), stringr::str_interp('${name}/${name}_cell_types.csv'))
  write.csv(table(obj$pred.type), stringr::str_interp('${name}/${name}_celltypes.csv'))
  
  Schwann <- subset(obj, subset =  pred.type == 'Schwann')
  Schwann <- FindVariableFeatures(Schwann)
  Schwann <- ScaleData(Schwann)
  Schwann <- RunPCA(Schwann, npcs = 30)
  Schwann <- RunUMAP(Schwann, reduction = 'pca', dims = 1:30)
  Schwann <- FindNeighbors(Schwann, reduction = 'pca', dims = 1:30)
  Schwann <- FindClusters(Schwann, resolution = 1)
  schwann_only <- DimPlot(Schwann, reduction = 'umap')
  ggsave(stringr::str_interp('${name}/${name}_schwann_only.png'), plot = schwann_only)
  
  schwann_features <- FeaturePlot(Schwann, features = c('CADM2', 'CDH19','FXYD1',
                                                        'LGI4','MAL','MPZ',
                                                        'PLP1','PMP2','PMP22','S100B'))
  ggsave(stringr::str_interp('${name}/${name}_schwann_features.png'), plot = schwann_features)
  #save object to unique name in the environment
  assign(name, obj)
}
paths <- c('scVS01.h5','scVS02.h5')

for (path in paths){
  name <- substr(path,1,6)
  filt.matrix <- Read10X_h5(path,use.names = T)
    names <- c(names, name)
  
  raw.matrix  <- Read10X_h5(stringr::str_interp("${name}_raw.h5"),use.names = T)
    #str(raw.matrix)
  #str(filt.matrix)
  
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  #soup.channel
  
  b <- CreateSeuratObject(counts = filt.matrix, assay = "RNA", project = stringr::str_interp("${name}"))
  b <- NormalizeData(b)
  b <- FindVariableFeatures(b)
  b <- ScaleData(b)
  b <- RunPCA(b)
  b <- RunUMAP(b, dims = 1:30)
  b <- FindNeighbors(b, dims = 1:30)
  b <- FindClusters(b, resolution = 0.5)
  meta <- b@meta.data
  umap <- b@reductions$umap@cell.embeddings
  rm(b)
  
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  soup.channel = autoEstCont(soup.channel)
  #ggsave(stringr::str_interp('${name}/${name}_soup_channel.png'), plot = autoEstCont(soup.channel))
  
  adjusted_counts <- adjustCounts(soup.channel)
  obj <- CreateSeuratObject(counts = adjusted_counts, assay = "RNA", project = stringr::str_interp("${name}"))
  obj$type ="NF2"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('${name}/${name}_QC_stats.png'), plot = QC_stats)
  #QC_stats
  ## Remove empty droplets and low quality cells
  obj <- subset(obj, nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 25)
  post_cutoff_QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('${name}/${name}_post_cutoff_QC_stats.png'), plot = post_cutoff_QC_stats)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  ## Scale
  obj <- ScaleData(obj, features = rownames(data))
  ## Run PCA
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = 1:30)
  ## Find and Remove Doublets
  # pK Identification (no ground-truth)
  sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  #ggsave(stringr::str_interp('${name}/${name}_pK_identify.png'), plot = bcmvn_obj)
  annotations <- Idents(obj)
  ## ex: annotations <- obj@meta.data$ClusteringResults
  
  homotypic.prop <- modelHomotypic(annotations)         
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  
  nExp_poi <- round(0.075*nrow(obj@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder
  obj <- doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) ## Unsure of most appropriate parameters, these are from the example on the doublet finder github
  # this column in meta.data have a different name for every sample and/or parameter setting, so I store the name of the column to refer to it later without hard coding
  
  doublet_col <- rev(names(obj@meta.data))[1] 
  doublets <- DimPlot(obj, reduction = 'umap', group.by =doublet_col)
  ggsave(stringr::str_interp('${name}/${name}_doublets.png'), plot = doublets)
  
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5)
  ## Visualize before analysis
  clusters <- DimPlot(obj, reduction = 'umap')
  
  ggsave(stringr::str_interp('${name}/${name}_clusters.png'), plot = clusters)
  
  markerlist <- FindAllMarkers(obj)
  write.csv(markerlist, file = stringr::str_interp("${name}/top_markerlist.csv"))
  
  ## Using scSorter
  topgenes <- head(VariableFeatures(obj), 2000)
  obj_exp = GetAssayData(obj)
  topgene_filter = rowSums(as.matrix(obj_exp)[topgenes, ]!=0) > ncol(obj_exp)*.1
  topgenes = topgenes[topgene_filter]
  ## load annotation file
  rm(anno)
  anno <- read.csv("sorted_sc_markers_DA.csv")
  anno
  picked_genes = unique(c(anno$Marker, topgenes))
  write.csv(picked_genes,file = stringr::str_interp("${name}/${name}_picked_genes.csv"))
  obj_exp = obj_exp[rownames(obj_exp) %in% picked_genes, ]
  ## run scSorter
  rts <- scSorter(obj_exp, anno)
  rts
  obj$pred.type <- rts$Pred_Type
  table(obj$pred.type)
  
  cell_types <- DimPlot(obj,  reduction = 'umap' , group.by = "pred.type")
  ggsave(stringr::str_interp('${name}/${name}_celltypes.png'), plot = cell_types)
  ## Tabulate Cell Types
  write.csv(table(Idents(obj)), stringr::str_interp('${name}/${name}_cell_types.csv'))
  write.csv(table(obj$pred.type), stringr::str_interp('${name}/${name}_celltypes.csv'))
  
  Schwann <- subset(obj, subset =  pred.type == 'Schwann')
  Schwann <- FindVariableFeatures(Schwann)
  Schwann <- ScaleData(Schwann)
  Schwann <- RunPCA(Schwann, npcs = 30)
  Schwann <- RunUMAP(Schwann, reduction = 'pca', dims = 1:30)
  Schwann <- FindNeighbors(Schwann, reduction = 'pca', dims = 1:30)
  Schwann <- FindClusters(Schwann, resolution = 1)
  schwann_only <- DimPlot(Schwann, reduction = 'umap')
  ggsave(stringr::str_interp('${name}/${name}_schwann_only.png'), plot = schwann_only)
  
  schwann_features <- FeaturePlot(Schwann, features = c('CADM2', 'CDH19','FXYD1',
                                                        'LGI4','MAL','MPZ',
                                                        'PLP1','PMP2','PMP22','S100B'))
  ggsave(stringr::str_interp('${name}/${name}_schwann_features.png'), plot = schwann_features)
  #save object to unique name in the environment
  assign(name, obj)
}
#does not need operator

paths <- c('scVS10.h5' ,'scVS11.h5', 'snVS10.h5' ,'snVS11.h5', 'snVS12.h5' ,'snVS13.h5','snVS14.h5','snVS15.h5')
for (path in paths){
  name <- substr(path,1,6)
  filt.matrix <- Read10X_h5(path,use.names = T)
#  filt.matrix <- filt.matrix$`Gene Expression`
  names <- c(names, name)
  
  raw.matrix  <- Read10X_h5(stringr::str_interp("${name}_raw.h5"),use.names = T)
 # raw.matrix <- raw.matrix$`Gene Expression`
  #str(raw.matrix)
  #str(filt.matrix)
  
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  #soup.channel
  
  b <- CreateSeuratObject(counts = filt.matrix, assay = "RNA", project = stringr::str_interp("${name}"))
  b <- NormalizeData(b)
  b <- FindVariableFeatures(b)
  b <- ScaleData(b)
  b <- RunPCA(b)
  b <- RunUMAP(b, dims = 1:30)
  b <- FindNeighbors(b, dims = 1:30)
  b <- FindClusters(b, resolution = 0.5)
  meta <- b@meta.data
  umap <- b@reductions$umap@cell.embeddings
  rm(b)
  
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  soup.channel = autoEstCont(soup.channel)
  #ggsave(stringr::str_interp('${name}/${name}_soup_channel.png'), plot = autoEstCont(soup.channel))
  
  adjusted_counts <- adjustCounts(soup.channel)
  obj <- CreateSeuratObject(counts = adjusted_counts, assay = "RNA", project = stringr::str_interp("${name}"))
  obj$type ="Sporadic"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('${name}/${name}_QC_stats.png'), plot = QC_stats)
  #QC_stats
  ## Remove empty droplets and low quality cells
  obj <- subset(obj, nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 25)
  post_cutoff_QC_stats <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(stringr::str_interp('${name}/${name}_post_cutoff_QC_stats.png'), plot = post_cutoff_QC_stats)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  ## Scale
  obj <- ScaleData(obj, features = rownames(data))
  ## Run PCA
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = 1:30)
  ## Find and Remove Doublets
  # pK Identification (no ground-truth)
  sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  #ggsave(stringr::str_interp('${name}/${name}_pK_identify.png'), plot = bcmvn_obj)
  annotations <- Idents(obj)
  ## ex: annotations <- obj@meta.data$ClusteringResults
  
  homotypic.prop <- modelHomotypic(annotations)         
  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  
  nExp_poi <- round(0.075*nrow(obj@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder
  obj <- doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE) ## Unsure of most appropriate parameters, these are from the example on the doublet finder github
  # this column in meta.data have a different name for every sample and/or parameter setting, so I store the name of the column to refer to it later without hard coding
  
  doublet_col <- rev(names(obj@meta.data))[1] 
  doublets <- DimPlot(obj, reduction = 'umap', group.by =doublet_col)
  ggsave(stringr::str_interp('${name}/${name}_doublets.png'), plot = doublets)
  
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5)
  ## Visualize before analysis
  clusters <- DimPlot(obj, reduction = 'umap')
  
  ggsave(stringr::str_interp('${name}/${name}_clusters.png'), plot = clusters)
  
  markerlist <- FindAllMarkers(obj)
  write.csv(markerlist, file = stringr::str_interp("${name}/top_markerlist.csv"))
  
  ## Using scSorter
  topgenes <- head(VariableFeatures(obj), 2000)
  obj_exp = GetAssayData(obj)
  topgene_filter = rowSums(as.matrix(obj_exp)[topgenes, ]!=0) > ncol(obj_exp)*.1
  topgenes = topgenes[topgene_filter]
  ## load annotation file
  rm(anno)
  anno <- read.csv("sorted_sc_markers_DA.csv")
  anno
  picked_genes = unique(c(anno$Marker, topgenes))
  write.csv(picked_genes,file = stringr::str_interp("${name}/${name}_picked_genes.csv"))
  obj_exp = obj_exp[rownames(obj_exp) %in% picked_genes, ]
  ## run scSorter
  rts <- scSorter(obj_exp, anno)
  rts
  obj$pred.type <- rts$Pred_Type
  table(obj$pred.type)
  
  cell_types <- DimPlot(obj,  reduction = 'umap' , group.by = "pred.type")
  ggsave(stringr::str_interp('${name}/${name}_celltypes.png'), plot = cell_types)
  ## Tabulate Cell Types
  write.csv(table(Idents(obj)), stringr::str_interp('${name}/${name}_cell_types.csv'))
  write.csv(table(obj$pred.type), stringr::str_interp('${name}/${name}_celltypes.csv'))
  
  Schwann <- subset(obj, subset =  pred.type == 'Schwann')
  Schwann <- FindVariableFeatures(Schwann)
  Schwann <- ScaleData(Schwann)
  Schwann <- RunPCA(Schwann, npcs = 30)
  Schwann <- RunUMAP(Schwann, reduction = 'pca', dims = 1:30)
  Schwann <- FindNeighbors(Schwann, reduction = 'pca', dims = 1:30)
  Schwann <- FindClusters(Schwann, resolution = 1)
  schwann_only <- DimPlot(Schwann, reduction = 'umap')
  ggsave(stringr::str_interp('${name}/${name}_schwann_only.png'), plot = schwann_only)
  
  schwann_features <- FeaturePlot(Schwann, features = c('CADM2', 'CDH19','FXYD1',
                                                        'LGI4','MAL','MPZ',
                                                        'PLP1','PMP2','PMP22','S100B'))
  ggsave(stringr::str_interp('${name}/${name}_schwann_features.png'), plot = schwann_features)
  #save object to unique name in the environment
  assign(name, obj)
}

# integrate ---------------------------------------------------------------
names
rm(clusters, doublets, bcmvn_obj, post_cutoff_QC_stats, QC_stats, sweep.res.list_obj, sweep.stats_obj)
dir.create('integrated/')

obj.list <- mget(names)

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated, features = rownames(data))
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = 'pca', dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = 'pca', dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)
p1 <- DimPlot(integrated, reduction = 'umap')
ggsave('umap_default_clusters.png', plot = p1)
Idents(integrated) <- 'pred.type'
# loop for DEGs for all cell-types versus integrated ---------
library(Seurat)

cell_types <- c("Schwann", "Macrophage", "Endothelial", "Neutrophil", "RBCs",
                "T-cells", "Proliferating Fibroblast", "Unknown", "Pericytes", "NK Cells")

# Define the Seurat object if it is not already defined
integrated <- Read10X_h5("integrated/integrated.h5")

# Loop through each cell type
for (cell_type in cell_types) {
  # Find DEGs for the current cell type

DEGs <- FindMarkers(integrated, ident.1 = Schwann)
  
  # Check if DEGs were found
  if (nrow(DEGs) > 0) {
    # Append DEGs to CSV file in a new sheet
    sheet_name <- paste0("DEGs_", cell_type)
    write.table(DEGs, file = paste0(file_path, ".", cell_type, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE, append = FALSE)
  } else {
    # Print a message if no DEGs were found
    cat("No DEGs found for", cell_type, "\n")
  }
}



p2 <- DimPlot(integrated, reduction = 'umap', group.by = 'pred.type')
ggsave('integrated/umap_clusters_scSorter.png', plot = p1+p2)
p1 + p2
ggsave('integrated/integrated_umap.png', plot = p1+p2)
table(Idents(integrated))

write.csv(table(Idents(integrated)), 'integrated/Idents_cell_types.csv')

Schwann <- subset(integrated, subset =  pred.type == 'Schwann')
DefaultAssay(Schwann) <- 'RNA'
write.csv(table(Idents(integrated),  integrated$orig.ident), 'integrated/Idents_orig.csv')
write.csv(table(Idents(integrated), integrated$type, integrated$orig.ident), 'integrated/Idents_orig_type.csv')
# Secondary Analysis (reclustering)
Schwann <- FindVariableFeatures(Schwann)
Schwann <- ScaleData(Schwann, features = row.names(Schwann))
Schwann <- RunPCA(Schwann, npcs = 30)
Schwann <- RunUMAP(Schwann, reduction = 'pca', dims = 1:30)
Schwann <- FindNeighbors(Schwann, reduction = 'pca', dims = 1:30)
Schwann <- FindClusters(Schwann, resolution = 0.1)
write.csv(FindAllMarkers(Schwann), 'schwann_markers.csv')

s1 <- DimPlot(Schwann, reduction = 'umap')
ggsave('integrated/umap_Schwann_only_clusters.png', plot = s1)
s2 <- DimPlot(Schwann, reduction = 'umap', split.by = 'type')
ggsave('integrated/umap_Schwann_only_by_ident.png', plot = s2)
s3 <-  DimPlot(Schwann, reduction = 'umap', split.by = 'orig.ident', ncol =5)
ggsave('integrated/umap_Schwann_orig_ident.png', plot = s3)
p3 <- DimPlot(integrated, reduction = 'umap', split.by = 'orig.ident',ncol = 5)
ggsave('integrated/umap_clusters_scSorter.png', plot = p3)
p4 <- DimPlot(integrated, reduction = 'umap', split.by = 'type')
ggsave('integrated/umap_NF2_vs_Sporadic.png', plot = p4)

rm(Sporadic,NF2)

Sporadic <- subset(integrated, subset = type == 'Sporadic')

s_umap <- DimPlot(Sporadic,split.by = 'orig.ident', ncol=3)

ggsave('integrated/umap_clusters_Sporadic.png', plot = s_umap)
s_features <- FeaturePlot(Sporadic, features = c("VEGFA","TEAD1","YAP1"))
ggsave('integrated/umap_features_Sporadic.png', plot = s_features)

so_features <- FeaturePlot(Sporadic, split.by = 'orig.ident', features = c("VEGFA","TEAD1","YAP1"))
ggsave('integrated/umap_features_Sporadic_split.png', plot = so_features)

NF2 <- subset(integrated, subset = type == 'NF2')

n_umap <- DimPlot(NF2,split.by = 'orig.ident', ncol=4)
ggsave('integrated/umap_clusters_NF2.png', plot = n_umap)

n_features <- FeaturePlot(NF2, features = c("VEGFA","TEAD1","YAP1"))
ggsave('integrated/umap_features_NF2.png', plot = n_features)

no_features <- FeaturePlot(NF2, split.by = 'orig.ident', features = c("VEGFA","TEAD1","YAP1"))
ggsave('integrated/umap_features_NF2_split.png', plot = no_features)
ggsave('integrated/umap_clusters_NF2_Sporadic.png', plot = n_umap + s_umap)




Macrophage <- subset(integrated, subset =  pred.type == 'Macrophage')
DefaultAssay(Macrophage) <- 'RNA'
# Secondary Analysis (reclustering)
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, features = row.names(Macrophage))
Macrophage <- RunPCA(Macrophage, npcs = 30)
Macrophage <- RunUMAP(Macrophage, reduction = 'pca', dims = 1:30)
Macrophage <- FindNeighbors(Macrophage, reduction = 'pca', dims = 1:30)
Macrophage <- FindClusters(Macrophage, resolution = 0.1)
write.csv(FindAllMarkers(Macrophage), 'macrophage_markers.csv')

m1 <- DimPlot(Macrophage, reduction = 'umap')
ggsave('integrated/umap_Macrophage_only_clusters.png', plot = m1)
m2 <-  DimPlot(Macrophage, reduction = 'umap', split.by = 'type')
ggsave('integrated/umap_Macrophage_type.png', plot = m2)
m3 <-  DimPlot(Macrophage, reduction = 'umap', split.by = 'orig.ident', ncol =5)
ggsave('integrated/umap_Macrophage_orig_ident.png', plot = m3)

save.image("Renvironment.Rdata")
saveRDS(Sporadic,file = "Robject_Sporadic.rds")
saveRDS(NF2,file = "Robject_NF2.rds")
saveRDS(Macrophage,file = "Robject_Macrophage.rds")
saveRDS(Schwann,file = "Robject_Schwann.rds")
saveRDS(integrated,file = "Robject_integrated.rds")
for (name in names) {
  obj <- get(name)  
  saveRDS(obj, file = paste0(name, ".rds"))  # Save the object to a file
}
