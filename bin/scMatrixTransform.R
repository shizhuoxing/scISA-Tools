rm(list=ls())
library(dplyr)
library(Seurat)
library(Matrix)

args<-commandArgs(TRUE)
data_dir<-args[1]
out_dir<-args[2]
prefix<-args[3]
trtype<-args[4]

if(trtype=="text2sparse"){library(DropletUtils)}

if(trtype!="text2sparse"){metadata <- Read10X(data.dir =data_dir)}else{metadata <- read.csv(data_dir, row.names =1, sep=",")}
if(trtype=="text2sparse"){data <- CreateSeuratObject(counts = metadata, project = prefix)}
string<-paste(out_dir,prefix,".matrix.csv",sep="")
if(trtype!="text2sparse"){write.table(metadata, file = string, sep =",", row.names =TRUE, col.names =TRUE, quote =TRUE)}else{write10xCounts(x =data@assays$RNA@counts, path = out_dir)}
