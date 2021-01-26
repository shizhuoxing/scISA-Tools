rm(list=ls())
library(dplyr)
library(Seurat)
library(Matrix)
library(argparse)
library(stringr)

args<-commandArgs(TRUE)
data_dir<-args[1]
sample_name<-args[2]
outdir<-args[3]
type<-args[4]

if(type!="csv"){ml6tpen.metadata <- as.matrix(Read10X(data.dir =data_dir))}else{ml6tpen.metadata <- read.csv(data_dir, row.names =1, sep=",")}

ml6tpen <- CreateSeuratObject(counts = ml6tpen.metadata, min.cells = 5, min.features = 100, project = sample_name)
pbmc <- NormalizeData(ml6tpen, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

mat<-Embeddings(pbmc,"tsne")
mat<-as.data.frame(mat)
mat$Barcode<-rownames(mat)
mat<-mat[,c(3,1,2)]
string<-paste(outdir,sample_name,".tsne_projection.csv",sep="")
write.table(mat,file=string,sep=",",quote=FALSE,row.names=FALSE)

mat<-Embeddings(pbmc,"umap")
mat<-as.data.frame(mat)
mat$Barcode<-rownames(mat)
mat<-mat[,c(3,1,2)]
string<-paste(outdir,sample_name,".umap_projection.csv",sep="")
write.table(mat,file=string,sep=",",quote=FALSE,row.names=FALSE)

mat<-as.data.frame(pbmc@active.ident)
string<-paste(outdir,sample_name,".cellBC_cluster.csv",sep="")
write.table(mat,file=string,sep=",",quote=FALSE,row.names=TRUE, col.names =FALSE)
