#cleanup
#rm(meta, Data)
#Clustering Via SNN
Object <- FindNeighbors(Object, dims = 1:15)
for(r in seq(0, 2, by = 0.2)){
  Object <- FindClusters(Object, resolution = r, algorithm = 1)
}

#Application of cluste tree to find optimal number of clusters
clustree(Object)
Object <- FindClusters(Object, resolution = 0.2, algorithm = 1)

#NMI and ARI for clustering and metadata
Obj_Cluster <- clustComp(Idents(Object),Object@meta.data$cluster_label)
Obj_Class <- clustComp(Idents(Object),Object@meta.data$class_label)

#Running UMAP and TSNE
Object <- RunUMAP(Object, dims = 1:15)
Object <- RunTSNE(Object, dims = 1:15)

#Plotting the clustering
DimPlot(Object, reduction = "umap", label = TRUE)
DimPlot(Object, reduction = "tsne", label = TRUE)

#Extraction of cluster markers
Object.markers <- FindAllMarkers(Object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Grouping Clusters
Object.markers %>%
  group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(Object, features = top10)

#Top 10 genes in each cluster
Object.markers %>%
  group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(Object, features = top10$gene) + NoLegend()

#Extraction of gene list
cluster1.markers<- FindMarkers(Object, ident.1 = 1)
n.gt.0 <- apply(cluster1.markers, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(cluster1.markers)[which(n.gt.0/ncol(cluster1.markers) >= 0.75)]
all.genes <- rownames(cluster1.markers)
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

#Application of GO
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")        
GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)