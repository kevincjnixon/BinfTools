## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=8,
  fig.height=8,
  fig.align = "center"
)

## ----setup--------------------------------------------------------------------
#devtools::install_github("kevincjnixon/BinfTools")
library(BinfTools)

## ----data---------------------------------------------------------------------
###Load data generated from DESeq2:
head(norm_counts) #Normalized gene counts (from counts(dds, normalized=T))
head(res) #Results object (from as.data.frame(results(dds, contrast=c("Condition","WT","KO"))))
print(cond) #Character vector of conditions (from as.character(dds$Condition))
print(geneSet) #List of gene sets related to rhodopsin signaling imported from qusage::read.gmt("rhodopsin.gmt")

## ----converting---------------------------------------------------------------
### Converting Limma results:
#res<-limma::topTable(fit, number=Inf)
#res<-fromLimma(res)
###Converting EdgeR results:
#res<-edgeR::topTags(de, n=Inf)
#res<-fromEdgeR(res)

## ----GSEA, warning=FALSE------------------------------------------------------
GenerateGSEA(res, filename="GSEA.rnk", bystat=T, byFC=F)

## ----MA_Plot, fig.cap="MA Plot", message=FALSE--------------------------------
#Pull top significant genes to label and colour
genes<-rownames(res[order(res$padj),])[1:5]
#Make an MA plot
MA_Plot(res, title="MA Plot", p=0.05, FC=log(1.25,2), lab=genes, col=genes)

## ----volcano, fig.cap="Volcano Plot diplaying upregulated (red) and downregulated (blue) genes", message=FALSE----
volcanoPlot(res, title="Volcano Plot", p=0.05, FC=log(1.25,2), lab=genes, col=genes)

## ----GO, warning=FALSE, fig.cap=c("Top 10 Enriched GO terms - minimum 10 genes/term","Top 10 Significant GO terms - maximum 500 genes/term")----
#Get downregulated gene names
genes<-rownames(subset(res, padj<0.05 & log2FoldChange < log(1.25, 2)))
#create output directory
dir.create("GO")
#Run GO analysis for Biological Process and output to folder "GO/GO_analysis*"
GO_GEM(genes, species="dmelanogaster", bg=rownames(res), source="GO:BP", prefix="GO/GO_analysis")

## ----heatmap, fig.cap="Heatmap of DEGs"---------------------------------------
#Get the names of all differentially expressed genes
DEGs<-rownames(subset(res, padj<0.05 & abs(log2FoldChange)>log(1.25,2)))
#Pull top significant genes to label and colour
top_genes<-rownames(res[order(res$padj),])[1:5]
#Generate a heatmap
zheat(genes=DEGs, counts=norm_counts, conditions=cond, con="WT", title="DEGs", labgenes = top_genes)

## ----count_plot, fig.cap="Violin plot comparing gene expression of upregulated genes between conditions"----
#Check all upregulated genes
genes<-rownames(subset(res, padj<0.05 & log2FoldChange > log(1.25,2)))
#Compare gene expression between the two conditions using z-score normalized counts:
count_plot(counts=norm_counts, scaling="zscore", genes=genes, condition=cond, title="Upregulated Genes", compare=list(c("WT","KO")))

## ----gsva, fig.cap="Violin plot comparing normalized enrichment scores of pathways involved in rhodopsin signaling between conditions"----
#Import the gene sets associated with pathways of an enrichment map cluster - In this case rhodopsin signaling
#library(qusage)
#geneSet<-read.gmt("rhodopsin.gmt")

#Now run the gsva using z-score normalized gene expression and make a plot:
gsva_plot(counts=t(scale(t(norm_counts))), geneset=geneSet, method="ssgsea", condition=cond, title="Rhodopsin-Mediated Signaling", compare=list(c("WT","KO")))


