library(Seurat)
library(sceasy)

#convert H5AD
sceasy::convertFormat("path/to7object.h5ad",from="anndata", to="seurat",outFile='output/object.rds')

#read in and ensure filtered
limb1<-readRDS("output/object.rds")
limb1 <- subset(x = limb1, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent_mito < 0.1)

#import concatenated scrublet outputs#
doublets_limb <- read.table("path/to/scublet_output.tsv",header = F,row.names = 1)
colnames(doublets_limb) <- c("Doublet_score","Is_doublet")
limb1$Doublet_score<- doublets_limb$Doublet_score
limb1$Is_doublet<- doublets_limb$Is_doublet

#plot qc
VlnPlot(limb1, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb", "Doublet_score"),ncol = 5)

#run UMAP and overcluster for doublet detection
limb1 <- RunUMAP(object = limb1, dims = 1:50,min.dist = 0.5,n.neighbors = 50, verbose = T)
limb1 <- FindNeighbors(object = limb1, dims = 1:20, verbose = T)
limb1 <- FindClusters(object = limb1, resolution = 10, verbose = T)
table(limb1@active.ident, limb1$Is_doublet) 

#remove any clusters >20% consisting of doublets
limb1<-subset(limb1, idents="189", invert=T)

#remove individual doublet cells
limb1<-SetIdent(limb1, value = "Is_doublet")
limb1<-subset(limb1, idents = "False")
