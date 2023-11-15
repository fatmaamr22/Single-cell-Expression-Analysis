#Reading .csv Data and attaching metadata
Data <- read_csv("../Datasets/Alan_mouse_LGN/mouse_LGN_2021_exon-matrix.csv", col_names = TRUE)
Data <- as.data.frame(Data)
rownames(Data) <- Data[, 1]
Data <- Data[, -1]
meta <- as.data.frame(read.csv("../Datasets/Alan_mouse_LGN/mouse_LGN_2021_metadata.csv",row.names = 2,header = TRUE))
meta <- meta[,-1]
meta <- meta[rownames(meta) %in% colnames(Data),]

#Creating seurat object
Object <- CreateSeuratObject(counts = Data,names.delim = "-",meta.data = meta, names.field = 3)

#VLN and FScatter plots
VlnPlot(Object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(Object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Normalization and variable features extraction
Object <- NormalizeData(Object)
Object <- FindVariableFeatures(Object, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Object), 10)

#Variable features plot
plot <- VariableFeaturePlot(Object)
LabelPoints(plot = plot, points = top10, repel = TRUE, xnudge =  0 , ynudge = 0)

#PCA
all.genes <- rownames(Object)
Object <- ScaleData(Object, features = all.genes)
Object <- RunPCA(Object, features = VariableFeatures(object = Object))

#PCA graph
print(Object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Object, dims = 1:2, reduction = "pca")
DimHeatmap(Object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Object, dims = 1:15, cells = 500, balanced = TRUE)

#preforming Jackstraw analysis
Object <- JackStraw(Object, num.replicate = 100)
Object <- ScoreJackStraw(Object, dims = 1:20)

#Plotting Jackstraw graph and elbow plots
JackStrawPlot(Object, dims = 1:20)
ElbowPlot(Object)