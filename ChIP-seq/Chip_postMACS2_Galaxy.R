setwd("~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/data_analysis/R_analysis/VS/results/Chip")

#Download MACS2 callpeaks tabular form from Galaxy
macsPeaks <- "MACS2_callpeaks.xls"
macsPeaks_DF <- read.delim(macsPeaks, comment.char = "#")
macsPeaks_DF[1:2, ]

#Get GRanges to get information of the peaks from the pulldown 
library(GenomicRanges)
macsPeaks_GR <- GRanges(seqnames = macsPeaks_DF[, "chr"], IRanges(macsPeaks_DF[,"start"], macsPeaks_DF[, "end"]))
macsPeaks_GR
seqnames(macsPeaks_GR)
mcols(macsPeaks_GR) <- macsPeaks_DF[, c("abs_summit", "fold_enrichment")]
macsPeaks_GR

#annotate the peaks 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomeInfoDb)
library(ChIPseeker)

#For 500 bp before and after the center of the pulldown 
{
peakAnno_500 <- annotatePeak(macsPeaks_GR, tssRegion = c(-500, 500), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                         annoDb = "org.Hs.eg.db")
peakAnno_500
annotatedPeaksGR <- as.GRanges(peakAnno_500)
annotatedPeaksDF <- as.data.frame(peakAnno_500)
annotatedPeaksDF[1:2, ]

#Plot location of genes relative to TSS
plotAnnoBar(peakAnno_500)
plotDistToTSS(peakAnno_500)

#Identify annotation peaks in the promoter regions
annotatedPeaksGR_TSS_500 <- annotatedPeaksGR[annotatedPeaksGR$annotation == "Promoter",
                                         ]
genesWithPeakInTSS_500 <- unique(annotatedPeaksGR_TSS_500$geneId)
genesWithPeakInTSS_500[1:2]

#Get gene names from hg38 database 
allGeneGR <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)
allGeneGR[1:2, ]
allGeneIDs <- allGeneGR$gene_id

#GSEA
library(clusterProfiler)
GO_result_500 <- enrichGO(gene = genesWithPeakInTSS_500, universe = allGeneIDs, OrgDb = org.Hs.eg.db,
                      ont = "BP")
GO_result_df_500 <- data.frame(GO_result_500)
GO_result_df[1:10, ]
#network plot
library(enrichplot)
GO_result_plot_500 <- pairwise_termsim(GO_result_500)
emapplot(GO_result_plot_500, showCategory = 20)

#Using another database to get the GO terms (experimental)
library(msigdbr)
msig_collections <- msigdbr_collections()
head(msig_collections, 25)
msig_t2g <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = NULL)
msig_t2g <- msig_t2g[, colnames(msig_t2g) %in% c("gs_name", "entrez_gene")]
msig_t2g[1:3, ]
hallmark <- enricher(gene = genesWithPeakInTSS, universe = allGeneIDs, TERM2GENE = msig_t2g)
hallmark_df <- data.frame(hallmark)
hallmark_df[1:5, ]
hallmark_plot <- pairwise_termsim(hallmark)
emapplot(hallmark_plot, showCategory = 20)

}

#For 100 bp before and after the center of the pulldown 
{
peakAnno_100 <- annotatePeak(macsPeaks_GR, tssRegion = c(-100, 100), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                         annoDb = "org.Hs.eg.db")
peakAnno_100
annotatedPeaksGR <- as.GRanges(peakAnno_100)
annotatedPeaksDF <- as.data.frame(peakAnno_100)
annotatedPeaksDF[1:2, ]
plotAnnoBar(peakAnno_100)
plotDistToTSS(peakAnno_100)
annotatedPeaksGR_TSS_100 <- annotatedPeaksGR[annotatedPeaksGR$annotation == "Promoter",
]
genesWithPeakInTSS_100 <- unique(annotatedPeaksGR_TSS_100$geneId)
genesWithPeakInTSS_100[1:2]
allGeneGR <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)
allGeneGR[1:2, ]
allGeneIDs <- allGeneGR$gene_id
library(clusterProfiler)
GO_result_100 <- enrichGO(gene = genesWithPeakInTSS_100, universe = allGeneIDs, OrgDb = org.Hs.eg.db,
                      ont = "BP")
GO_result_100_df <- data.frame(GO_result_100)
GO_result_100_df[1:10, ]
}

#For 750 bp before and after the center of the pulldown 
{

peakAnno_750 <- annotatePeak(macsPeaks_GR, tssRegion = c(-750, 750), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                             annoDb = "org.Hs.eg.db")
peakAnno_750
annotatedPeaksGR <- as.GRanges(peakAnno_750)
annotatedPeaksDF <- as.data.frame(peakAnno_750)
annotatedPeaksDF[1:2, ]
plotAnnoBar(peakAnno_750)
plotDistToTSS(peakAnno_750)
annotatedPeaksGR_TSS_750 <- annotatedPeaksGR[annotatedPeaksGR$annotation == "Promoter",
]
genesWithPeakInTSS_750 <- unique(annotatedPeaksGR_TSS_750$geneId)
genesWithPeakInTSS_750[1:2]
allGeneGR <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)
allGeneGR[1:2, ]
allGeneIDs <- allGeneGR$gene_id
library(clusterProfiler)
GO_result_750 <- enrichGO(gene = genesWithPeakInTSS_750, universe = allGeneIDs, OrgDb = org.Hs.eg.db,
                          ont = "BP")
GO_result_750_df <- data.frame(GO_result_100)
GO_result_750_df[1:10, ]
}

#For 1000 bp before and after the center of the pulldown 
{
peakAnno_1000 <- annotatePeak(macsPeaks_GR, tssRegion = c(-100, 100), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                              annoDb = "org.Hs.eg.db")
peakAnno_1000
annotatedPeaksGR <- as.GRanges(peakAnno_1000)
annotatedPeaksDF <- as.data.frame(peakAnno_1000)
annotatedPeaksDF[1:2, ]
plotAnnoBar(peakAnno_1000)
plotDistToTSS(peakAnno_1000)

annotatedPeaksGR_TSS_1000 <- annotatedPeaksGR[annotatedPeaksGR$annotation == "Distal Intergenic", ]
genesWithPeakInTSS_1000 <- unique(annotatedPeaksGR_TSS_1000$geneId)
genesWithPeakInTSS_1000[1:2]

allGeneGR <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)
allGeneGR[1:2, ]
allGeneIDs <- allGeneGR$gene_id

library(clusterProfiler)

# KEGG enrichment analysis
KEGG_result_1000 <- enrichKEGG(gene = genesWithPeakInTSS_1000, organism = "hsa", pvalueCutoff = 0.4, qvalueCutoff = 0.9)
KEGG_result_1000_df <- data.frame(KEGG_result_1000)
KEGG_result_1000_df[1:25, ]

# Plot KEGG enrichment results
pdf("KEGGTerms_25_DI.pdf", height = 20)
p1 <- dotplot(KEGG_result_1000, showCategory = 25)
print(p1)
dev.off()
entrez_genes <- mapIds(org.Hs.eg.db, genesWithPeakInTSS_1000, 'ENTREZID', 'SYMBOL')
# GSEA for KEGG
KEGG_result_1000_gsea <- gseKEGG(geneList = genesWithPeakInTSS_1000, organism = "hsa", nPerm = 1000)
KEGG_result_1000_gsea_df <- data.frame(KEGG_result_1000_gsea)
KEGG_result_1000_gsea_df[1:10, ]

# Plot GSEA results
pdf("KEGGTerms_GSEA_20.pdf")
p2 <- dotplot(KEGG_result_1000_gsea, showCategory = 20)
print(p2)
dev.off()
}
{
  peakAnno_1000 <- annotatePeak(macsPeaks_GR, tssRegion = c(-1000, 1000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                annoDb = "org.Hs.eg.db")
  peakAnno_1000
  annotatedPeaksGR <- as.GRanges(peakAnno_1000)
  annotatedPeaksDF <- as.data.frame(peakAnno_1000)
  annotatedPeaksDF[1:2, ]
  plotAnnoBar(peakAnno_1000)
  plotDistToTSS(peakAnno_1000)
  annotatedPeaksGR_TSS_1000 <- annotatedPeaksGR[annotatedPeaksGR$annotation == "Distal Intergenic",
  ]
  genesWithPeakInTSS_1000 <- unique(annotatedPeaksGR_TSS_1000$geneId)
  genesWithPeakInTSS_1000[1:2]
  allGeneGR <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)
  allGeneGR[1:2, ]
  allGeneIDs <- allGeneGR$gene_id
  library(clusterProfiler)
  GO_result_1000 <- enrichGO(gene = genesWithPeakInTSS_1000, universe = allGeneIDs, OrgDb = org.Hs.eg.db,
                             ont = "BP")
  GO_result_1000_df <- data.frame(GO_result_1000)
  GO_result_1000_df[1:10, ]
  GO_result_plot_1000 <- pairwise_termsim(GO_result_1000)
  emapplot(GO_result_plot_1000, showCategory = 20)
  pdf("GOTerms_30_DI.pdf", height = 20)
  p1 <- dotplot(GO_result_1000, showCategory=30)
  print(p1)
  dev.off()
  GO_result_1000 <- gseGO(gene = genesWithPeakInTSS_1000, universe = allGeneIDs, OrgDb = org.Hs.eg.db,
                          ont = "BP")
  GO_result_1000_df <- data.frame(GO_result_1000)
  GO_result_1000_df[1:10, ]
  GO_result_plot_1000 <- pairwise_termsim(GO_result_1000)
  emapplot(GO_result_plot_1000, showCategory = 20)
  pdf("GOTerms_20.pdf")
  p1 <- dotplot(GO_result_1000, showCategory=30)
  print(p1)
  dev.off()
}
#For 5000 bp before and after the center of the pulldown 
{
peakAnno_5000 <- annotatePeak(macsPeaks_GR, tssRegion = c(-5000, 5000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                              annoDb = "org.Hs.eg.db")
peakAnno_5000
annotatedPeaksGR <- as.GRanges(peakAnno_5000)
annotatedPeaksDF <- as.data.frame(peakAnno_5000)
annotatedPeaksDF[1:2, ]
plotAnnoBar(peakAnno_5000)
plotDistToTSS(peakAnno_5000)
annotatedPeaksGR_TSS_5000 <- annotatedPeaksGR[annotatedPeaksGR$annotation == "Promoter (<=1kb)",
]
genesWithPeakInTSS_5000 <- unique(annotatedPeaksGR_TSS_5000$geneId)
genesWithPeakInTSS_5000[1:2]
allGeneGR <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only=FALSE)
allGeneGR[1:2, ]
allGeneIDs <- allGeneGR$gene_id
library(clusterProfiler)
GO_result_5000 <- gseKEGG(geneList = genesWithPeakInTSS_5000, organism = "hsa", nPerm = 1000)
GO_result_5000_df <- data.frame(GO_result_5000)
GO_result_5000_df[1:10, ]
GO_result_plot_5000 <- pairwise_termsim(GO_result_5000)
emapplot(GO_result_plot_5000, showCategory = 20)
pdf("GOTerms_20.pdf")
p1 <- dotplot(GO_result_5000, showCategory=30)
print(p1)
dev.off()
}