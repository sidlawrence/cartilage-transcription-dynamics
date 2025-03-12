####cellsignalanalysis preperation and plotting####
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
source("scripts/cellSignalAnalysis_cellSignalAnalysis.R")
limb1<-readRDS("cartilage_v2/cartilage_v2.rds")
limb1<-CellCycleScoring(limb1, s.features= cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
limb1 <-SetIdent(limb1, value = limb1$Phase)
limb1<- subset(limb1, idents = "G1")

#biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
df<- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)

#create raw counts matrix
data<- limb1@assays$RNA@counts


#create rownames file
temp_rownames<-""
na_counter <- 1
for(i in 1:nrow(data)) {
  if(rownames(data)[i] %in% df$external_gene_name) {
    print(rownames(data)[i])
    print(df$ensembl_gene_id[which(df$external_gene_name == rownames(data)[i])] )
    temp_rownames[i] <- df$ensembl_gene_id[which(df$external_gene_name == rownames(data)[i])] 
  } else {
    print(paste(rownames(data)[i], " not found", sep = ""))
    temp_rownames[i] <- paste("NA", na_counter, sep = "")
    na_counter <- na_counter + 1
  }
}

rownames(data)<- temp_rownames
data<- data[!grepl("NA", rownames(data)),]

#remove duplicated rows
data<- data[!duplicated(rownames(data)), ]

#write rowNames file
rownames<- rownames(data)
write.table(rownames, file='cart_rowNames.tsv', sep='\t', row.names = F, quote = FALSE)

#writing corresponding column names
meta<- limb1@meta.data
columns<- meta$finalannoV3
write.table(columns, file='cart_columnNames.tsv', sep='\t', row.names = F, quote = FALSE)

#write matrix itself
writeMM(obj = data, file = "cart.mtx")

####see txt file to run####

####plotting####
fit<- normaliseExposures(tgt = "path/to/cellsignalanalysis/fitexposures.tsv")
plotExposures(fit = fit, cluster_rows = T, cluster_columns=F, row_names_gp = gpar(fontsize = 6))
