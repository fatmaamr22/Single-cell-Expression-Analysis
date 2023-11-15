#cleanup
#rm(resultFisher, GenTable, GOdata, cluster1.markers, Object.markers, plot, n.gt.0, expressed.genes)
#Extracting GO terms
geneList <- VariableFeatures(object = Object)
GOfilter <- groupGO(geneList,keyType = 'SYMBOL',
                    org.Mm.eg.db,ont = 'BP',
                    level = 3,
                    readable =  TRUE)
GOfilter <- GOfilter[GOfilter@result$Count>0]

#setup for pooling the dataset
finalData <- matrix(nrow = nrow(GOfilter) ,
                    ncol = Object@assays$RNA@counts@Dim[2],
                    data=0,
                    dimnames = list(GOfilter$ID,
                                    unlist(Object@assays$RNA@data@Dimnames[2])))
genescount <- GOfilter$Count

#pooling the dataset
Sys.time()
for(row in 1:nrow(finalData)){
  geneList <- GOfilter$geneID[as.integer(row)]
  genes <- str_split(geneList,'/', genescount[as.integer(row)])
  val <- Object@assays$RNA@counts[unlist(genes),]
  if(!is.null(nrow(val))){
    val <- colSums(val, na.rm = TRUE)    #summed
    #val <- colMeans(val, na.rm = TRUE)   #mean
  }
  finalData[row, ] <- val
}
Sys.time()

#Creating pooled Seurat obj
Object_Pooled<- CreateSeuratObject(finalData)

#Normalization
Object_Pooled <- NormalizeData(Object_Pooled)

#Finding variable features
Object_Pooled <- FindVariableFeatures(Object_Pooled, selection.method = "vst", nfeatures = 2000)
top10P <- head(VariableFeatures(Object_Pooled), 10)

#PCA
all.features <- rownames(Object_Pooled)
Object_Pooled <- ScaleData(Object_Pooled, features = all.features)
Object_Pooled <- RunPCA(Object_Pooled, features = VariableFeatures(object = Object_Pooled))

#PCA visualization
print(Object_Pooled[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Object_Pooled, dims = 1:2, reduction = "pca")
DimHeatmap(Object_Pooled, dims = 1:15, cells = 500, balanced = TRUE)

#Jackstraw analysis on pooled data
Object_Pooled <- JackStraw(Object_Pooled, num.replicate = 100)
Object_Pooled <- ScoreJackStraw(Object_Pooled, dims = 1:20)

#Plotting Jackstraw analysis and elbow plot
JackStrawPlot(Object_Pooled, dims = 1:20)
ElbowPlot(Object_Pooled)

#Clustering Via SNN
Object_Pooled <- FindNeighbors(Object_Pooled, dims = 1:15)
for(r in seq(0, 2, by = 0.2)){
  Object_Pooled <- FindClusters(Object_Pooled, resolution = r, algorithm = 1)
}

#Application of cluster tree to find optimal number of clusters
clustree(Object_Pooled)
Object_Pooled <- FindClusters(Object_Pooled, resolution = 0.2, algorithm = 1)

#NMI and ARI for clustering and metadata
ObjPooled_Cluster <- clustComp(Idents(Object_Pooled),Object@meta.data$cluster_label)
ObjPooled_Class <- clustComp(Idents(Object_Pooled),Object@meta.data$class_label)

#Running UMAP and TSNE
Object_Pooled <- RunUMAP(Object_Pooled, dims = 1:15)
Object_Pooled <- RunTSNE(Object_Pooled, dims = 1:15)

#Plotting the clustering
DimPlot(Object_Pooled, reduction = "umap", label = TRUE)
DimPlot(Object_Pooled, reduction = "tsne", label = TRUE)

#extracting cluster markers
Object_Pooled.markers <- FindAllMarkers(Object_Pooled,
                                        only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.1)
Object_Pooled.markers %>%
  group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

top10P <- head(Object_Pooled.markers$gene, 10)

DoHeatmap(Object_Pooled, features = top10P) + NoLegend()
FeaturePlot(Object_Pooled, features = top10P)
