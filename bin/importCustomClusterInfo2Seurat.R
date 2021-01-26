rm(list=ls())
library(dplyr)
library(Seurat)

args<-commandArgs(TRUE)
data_dir<-args[1]
sample_name<-args[2]
celltype<-args[3]
type<-args[4]

if(type!="csv"){ml6tpen.data <- Read10X(data.dir =data_dir)}else{ml6tpen.data <- read.csv(data_dir, row.names =1, sep=",")}
ml6tpen <- CreateSeuratObject(counts = ml6tpen.data, project = sample_name)
ml6tpen <- NormalizeData(ml6tpen, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(ml6tpen)
ml6tpen <- ScaleData(ml6tpen, features = all.genes)
ml6tpen <- FindVariableFeatures(ml6tpen, selection.method = "vst", nfeatures = 2000)
ml6tpen <- RunPCA(ml6tpen, features = VariableFeatures(object = ml6tpen))
ml6tpen <- JackStraw(ml6tpen, num.replicate = 100)
ml6tpen <- ScoreJackStraw(ml6tpen, dims = 1:20)
ml6tpen <- FindNeighbors(ml6tpen, dims = 1:10)
ml6tpen <- FindClusters(ml6tpen, resolution = 0.5)

# import custom cluster infomation for replace seurat clustering
a<-read.csv(celltype,header=F)
b<-as.factor(a[[2]])
c<-as.factor(a[[1]])
names(b)<-c
ml6tpen@active.ident<-b

ml6tpen <- RunUMAP(ml6tpen, dims = 1:10)
pdf("UMAPPlot.pdf")
DimPlot(ml6tpen, reduction = "umap")
dev.off()

ml6tpen.markers <- FindAllMarkers(ml6tpen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ml6tpen.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- ml6tpen.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
pdf(paste0(sample_name,".MarkerGenePlot.pdf"))
DoHeatmap(ml6tpen, features = top10$gene) + NoLegend()
dev.off()
