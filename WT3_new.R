####Seurat V3

####WT3_New

####Loading Needed Packages####

library(dplyr)
library(Seurat)
library(ggplot2)
####Loading 10X Dataset#####

# To change column name in RDS file#

WT_3_new@meta.data
tail(WT_3_new@meta.data)

WT_3_new@meta.data$orig.ident <- as.character(WT_3_new@meta.data$orig.ident)
WT_3_new@meta.data$orig.ident[WT_3_new@meta.data$orig.ident == "WT_3_new3k"] <- "WT"
WT_3_new@meta.data$newType <- WT_3_new@meta.data$orig.ident
WT_3_new@meta.data$newType [WT_3_new@meta.data$newType == "WT_3_new3k"] <- "WT"


WT_3_new.data <- Read10X(data.dir ="D:/10x-WT3-V3/filtered_feature_bc_matrix/")

####Initializing Seurat Object####

WT_3_new <- CreateSeuratObject(counts=WT_3_new.data, project= "WT_3_new3k", min.cells = 20, min.features = 200)
WT_3_new
#min.cells= a gene is found in at least "10" cells
#min.features= a cell containes a minimum of "200" features

###Quality Control####

###QC Mitochondrial Gene Expression
memory.limit()
memory.limit(size=128000)
#If have broken cell, will see increased percent of mt

WT_3_new[["percent.mt"]] <- PercentageFeatureSet(WT_3_new,pattern = "^mt-")

#Checking for Mt gene expression and adding data to matrix of WT_3_new
  #head(WT_3_new@meta.data,5)

###Visualizing nFeauture, nCount, and percent MT

VlnPlot(WT_3_new, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
#Saving as 1000xproportional

#boxplot(WT_3_new@meta.data$nFeature_RNA)
#maybe use boxplot to filter out cells with nFeatures considered outliers???

###FeatureScatter to see percentMito vs nCount and nFeature vs nCount

plot1 <- FeatureScatter(WT_3_new,feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(WT_3_new,feature1 = "nCount_RNA", feature2= "nFeature_RNA")

CombinePlots(plots=list(plot1, plot2))

###Filtering out based on nFeature and percentMt####

WT_3_new <- subset(WT_3_new, subset =nFeature_RNA >200 & nFeature_RNA <6000 & percent.mt <10)

#Filter to keep cells with >200 genes and less than 3600 genes and a percentMT <5
#median(WT_3_new@meta.data$nFeature_RNA)

####Normalizing Filtered Data####
##Log Normalization- normalizes feature expression measurements for each cell by the total expression and multiplies ration by scale factor of 10,000 by default
##Normalized values stored in WT_3_new[["RNA"]]@data

WT_3_new <- NormalizeData(WT_3_new, normalization.method = "LogNormalize", scale.factor = 10000)

####Identification of highly variable features####

WT_3_new <- FindVariableFeatures(WT_3_new, selection.method = "vst", nfeatures = 2000)

#Will select the top 2000 most variable genes in the dataset for downstream analysis

#Will show top 10 most highly variable genes

top10 <- head(VariableFeatures(WT_3_new),10)

#Plotting variable features with and without lables

plot1 <- VariableFeaturePlot(WT_3_new)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

####Scaling data in prep for PCA####

#Shifts expression so that hte mean expression across cells is 0

#Scales the expression of each gene so that the variance across cells is 1- giving equal weight in downstream analysis so that highly expressed genes do not skew/dominate the data

#Results stored in WT_3_new[["RNA"]]@scale.data
#Scaling all genes so that if do heatmap of gene not in top 2000, will not be overexpressing
all.genes <- rownames(WT_3_new)
WT_3_new <- ScaleData(WT_3_new,features = all.genes)

####PCA- linear dimensional reduction####

WT_3_new <- RunPCA(WT_3_new,features = VariableFeatures(object = WT_3_new))

#####Heatmap####

DimHeatmap(WT_3_new, dims = 1:15, cells = 500, balanced = TRUE)
# dims= Number of PCs want to see
# cells= uses the 500 most extreme cells to visualize heatmaps- speeds up process

###Choosing Significant Principal Components for downstream analysis####

##Jackstraw Method
#Taking 1% of data and looking for significant PCs- repeats 100 times 

WT_3_new <- JackStraw(WT_3_new, num.replicate = 100)
WT_3_new <- ScoreJackStraw(WT_3_new, dims = 1:20)

JackStrawPlot(WT_3_new, dims = 1:20)

##ElbowPlot

ElbowPlot(WT_3_new)

#####Cell Clustering#####

WT_3_new <- FindNeighbors(WT_3_new, dims = 1:20)
WT_3_new <- FindClusters(WT_3_new, resolution = 0.6)
#Will try and form clusters based on cells that have similar expression patterns
#resolution: 0.4-1.2 for 3k typically optimal
#increasing resolution increases the number of clusters

####UMAP#####

## Downloading UMAP package reticulate::py_install(packages = 'umap-learn')
WT_3_new <- RunUMAP(WT_3_new, dims = 1:15)
DimPlot(WT_3_new,reduction = "umap",label = TRUE)
#Saved as x by 390 height!
WT_3_new <- RunTSNE(WT_3_new, dims = 1:20)
DimPlot(WT_3_new,reduction = "tsne",label = TRUE) 
DimPlot(WT_3_new,reduction = "tsne")
DimPlot(WT_3_new)

###Find corelation between clusters###
av.exp <- AverageExpression(WT_3_new)$RNA
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, c( "0","1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"))
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()




# Add cell number per cluster to cluster labels
cell.num <- table(Idents(WT_3_new))

ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)
UMAPPlot(object = WT_3_new, do.return = T, label.size = 5) +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")


####Naming Clusters#####

new.cluster.ids <- c("Cluster 0-Bipolar Cells (n=746)","Cluster 1-Bipolar Cells (n=619)","Cluster 2-RGC (n=371)","Cluster 3-Rods (n=349)","Cluster 4-MGCs (n=333)","Cluster 5-Horizontal Cells (n=316)","Cluster 6-Cones (n=225)","Cluster 7-Rods (n=172)","Cluster 8-Microglial Cells (n=169)","Cluster 9-RPC (n=154)","Cluster 10-Amacrine Cells (n=136)","Cluster 11-Cones (n=121)","Cluster 12-Bipolar Cells (n=100)","Cluster 13-RPE (n=81)","Cluster 14-Bipolar Cells (n=69)","Cluster 15-Oligodendrocytes (n=48)","Cluster 16-Jun,myc,fos genes (n=12)")
new.cluster.ids <- c("BP-1","BP-2","Undefined","RPE-1","AC-1","Cones-1","AC-2","Cones-2","Committed neuronal cellS","Rods","HC","RGC", "Cones-3", "Erythroid cells","RPE-2","Rods-2","MG-1", "AC-3","AC-4","Microglia-1","MG-2","RPC","BP-3","AC-5","Oligodendrocytes","Microglia-2","AC-6", "AC-7")	
new.cluster.ids <- c("Bipolar cells","Bipolar cells","Undefined","RPE","Amacrine cells","Cones","Amacrine cells","Cones","Committed neuronal cellS","Rods","HC","RGC", "Cones", "Erythroid cells","RPE","Rods","Muller glial cells","Amacrine cells","Amacrine cells","Microglia","Muller glial cells","RPC","Bipolar cells","Amacrine cells","Oligodendrocytes","Microglia","Amacrine cells", "Amacrine cells")		
names(new.cluster.ids) <- levels(WT_3_new)
WT_3_new <- RenameIdents(WT_3_new, new.cluster.ids)
DimPlot(WT_3_new, reduction = "umap", label = TRUE, pt.size = 1.0) + NoLegend()
DimPlot(WT_3_new,reduction = "tsne",label = TRUE)






#change cluster names
WT_3_new <- RenameIdents( object = WT_3_new, '15' = '9')


##All markers
WT_3_new.markers <- FindAllMarkers(WT_3_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WT_3_new.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
write.csv(WT_3_new.markers,file = "D:/10x-WT3-V3/New folder/min cell 20/WT- 0.6res tsne markers.csv")

###Saving as an object#####

saveRDS(WT_3_new, file = "D:/10x-WT3-V3/New folder/min cell 20/WT_3_new tsne.rds")


##HeatMap#####

top10<- WT_3_new.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC)
DoHeatmap(WT_3_new, features = top10$gene) + NoLegend()
write.csv(top100,file = "C:/Users/abisa/Desktop/10x-WT3-V3/New folder/min cell 50/top100marker.csv")
top10<- WT_3_new.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
write.csv(top10,file = "D:/10x-WT3-V3/New folder/min cell 20/top10 markers.csv")

#subset to remove undefined cluster###
WT_3_new2 <- subset(WT_3_new,idents = c("0", "1", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23") )

new.cluster.ids <- c("0","1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")	

WT_3_new2@meta.data$new.ident <- Idents(WT_3_new2)
WT_3_new2 <- SetIdent(WT_3_new2, value = WT_3_new2$seurat_clusters)
WT_3_new2@meta.data
names(new.cluster.ids) <- levels(WT_3_new2)
WT_3_new2 <- RenameIdents(WT_3_new2, new.cluster.ids)
DimPlot(WT_3_new2,reduction = "tsne",label = TRUE, pt.size=0.5)
WT_3_new2.markers <- FindAllMarkers(WT_3_new2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(WT_3_new2.markers,file = "D:/10x-WT3-V3/New folder/min cell 20/WT- 0.6res new2 tsne undefined removed markers.csv")
new.cluster.ids <- c('0' ="Bipolar cells-1", '1' = "Rods",'3' ="Retinal Pigment epithelial cells-1",'4' ="Amacrine cells-1",'5' ="Cones-1",'6' = "Amacrine cells-2",'7' ="Cones-2",'8' ="Horizontal cells",'9' ="Retinal Ganglion cells",'10' ="Bipolar cells-2",'11' = "Muller glial cells",'12' = "Erythroid cells",'13' ="Neurogenic progenitor cells",'14' ="Retinal Pigment epithelial cells-2",'15' = "Amacrine cells-3" ,'16' ="Retinal Progenitor cells",'17' = "Muller glial cells",'18' ="Microglia-1",'19' ="Bipolar cells-3",'20' ="Amacrine cells-4",'21' ="Oligodendrocytes",'22' ="Microglia-2",'23' = "Amacrine cells-5")		
names(new.cluster.ids) <- levels(WT_3_new2)
WT_3_new2 <- RenameIdents(WT_3_new2, new.cluster.ids)
DimPlot(WT_3_new2,reduction = "tsne",label = TRUE, pt.size=0.5)
#colnames(WT_3_new2@meta.data)[new.ident] <- "seurat_clusters"#
WT_3_new2@meta.data
WT_3_new2@meta.data$seurat_clusters <- Idents(WT_3_new2)
WT_3_new2@meta.data$seurat_clusters <- new.integrated$new.ident
WT_3_new2$new.ident <- NULL
WT_3_new2@meta.data$orig.ident[WT_3_new2@meta.data$orig.ident == "WT_3_new3k"] <- "WT"
WT_3_new2$seurat_clusters<- as.character(WT_3_new2$seurat_clusters)
WT_3_new2$orig.ident<- as.character(WT_3_new2$orig.ident)
saveRDS(WT_3_new2, file = "C:/Users/abisa/Downloads/shiny WT/WT_3_new2.rds")

top10<- WT_3_new2.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(WT_3_new2, features = top10$gene) + NoLegend()
write.csv(top10,file = "D:/10x-WT3-V3/New folder/min cell 50/top10a markers.csv")

top5<- WT_3_new2.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
DoHeatmap(WT_3_new2, features = top5$gene) + NoLegend()

#####Feature Plots#####

FeaturePlot(WT_3_new, features = c("si:ch211-251b21.1","dut","aurkb","cdk1","top2a","rrm2","rbbp4","psat1","stmn1a","insm1a","fabp7a","pcna","marcksl1b","marcksb","marcksl1a","marcksa","ccna2","ccnb1"))
FeaturePlot(WT_3_new, features = c("arl3l1", "arl3l2", "ncaldb", "sema7a", "unc119b", "clul1") ,min.cutoff ="q9", pt.size=0.5)
FeaturePlot(WT_3_new, features = c("c1qa","c1qb","c1qbp","grna") ,min.cutoff ="q9", pt.size=0.5)
FeaturePlot(WT_3_new, features = c("pde6c", "arr3a", "arr3b", "clul1", "gnat2") ,min.cutoff ="q9", pt.size=0.5)
FeaturePlot(WT_3_new, features = c("vsx1", "rs1a", "lrit1a", "fezf2", "prkca", "rdh10a"),min.cutoff ="q9", pt.size=0.5)
FeaturePlot(WT_3_new, features = c("ccnd1", "myca", "aurkb", "stmn1a") ,min.cutoff ="q9", pt.size=0.5)
VlnPlot(object = WT_3_new, features = "pde6c", split.by = "orig.ident", pt.size = 0.0)

features <- c ("kidins220a", "lhx9", "pfn1", "cd59","chgb", "dgkaa", "her4.3", "cahz", "apoc1", "slc6a9", "slc18a3a", "sagb", "rpe65a", "creg1", "pde6c", "rbpms2a", "aqp9a", "cngb1a", "marcksl1a", "clul1", "syt2a", "kera", "gad2", "dct", "fabp11b", "rs1a", "cabp2a")
features <- c ("lhx9", "ccr9a", "cd59", "chgb", "prkca", "cahz", "apoc1", "slc6a9", "slc18a3a", "sagb", "rpe65a", "creg1", "ptgdsb.2", "rbpms2a", "rprmb", "marcksl1a", "rom1a", "pde6c", "gad2", "pax10", "dct", "clul1", "cabp2a")
features <- c ("slc6a9", "gad2", "chata", "chgb", "calb1", "slc18a3a")
a <- DotPlot(WT_3_new2, features = features )
a  + coord_flip() +theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
 

DotPlot(WT_3_new, features = c("slc6a9", "gad2", "chata", "chgb", "calb1", "slc18a3a")) + RotatedAxis()
DotPlot(WT_3_new, features = c("pde6b", "gnat1", "cnga1", "rom1a")) 
DotPlot(WT_3_new, features = c("opn1sw1", "opn1sw2", "opn1lw1", "opn1lw2", "opn1mw1", "opn1mw2", "pde6b"))
FeaturePlot(WT_3_new, features = c("slc6a9", "gad2", "chata", "chgb", "calb1", "slc18a3a") ,min.cutoff ="q9", pt.size=0.5)   
DoHeatmap(subset(WT_3_new, downsample = 100), features = c("opn1sw1", "opn1sw2", "opn1lw1", "opn1lw2", "opn1mw1", "opn1mw2"), size = 3)


avg <- AverageExpression(WT_3_new, features = NULL, add.ident = NULL, return.seurat = TRUE, verbose = TRUE)
DoHeatmap(avg, features = c("aurkb", "stmn1a", "dut", "ccnd1"), group.by = "ident", size = 3, angle = 0, combine = TRUE, draw.lines = FALSE) + theme(axis.text.x = element_text(size = 7), axis.line = element_line(colour = "#ffffff")) + scale_fill_gradient2(low = '#1000ff', mid = "#ffffff", high = '#aa0101', space = "Lab", na.value = "#ffffff", midpoint = 0, guide = "colourbar", aesthetics = "fill") + coord_flip(expand = TRUE, clip = "on")
DoHeatmap(avg, genes.use = c("ctnnb1", "yap1", "lats2", "lats1", "aurkb", "myca" ))  + scale_fill_gradient2(low = '#1000ff', mid = "#ffffff", high = '#aa0101', space = "Lab",  midpoint = 0, guide = "colourbar", aesthetics = "fill")
DoHeatmap(avg, genes.use = c("gst1", "sod1", "prdx2", "mt1", "gfap", "glul" ),  col.low = "#0000FF",
          col.mid = "#FFFAFA", col.high = "#FF0000" ) + theme_classic() 
cd_genes1 <- c("aurkb", "stmn1a", "fabp7a", "ccnd1", "ctnnb1", "myca")
cd_genes1 <- c("foxo1a","pik3ca","pik3r3a")
DoHeatmap(avg, features = cd_genes1, draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
####TSNE######

WT_3_new <- RunTSNE(WT_3_new, dims = 1:20)
DimPlot(WT_3_new, reduction = "tsne", label = TRUE)# + ggtitle(label = "FIt-SNE")


WT_3_new.markers <- FindAllMarkers(WT_3_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WT_3_new.markers %>% group_by(cluster) %>% top_n(n=2, wt=avg_logFC)
write.csv(WT_3_new.markers,file = "C:/Users/eshih/Desktop/WT3/WT3_New_TSNE_Cluster.csv")

subset9 <-subset(WT_3_new,idents = c("Rods"))
subset16 <-subset(WT_3_new,idents = c("MG-1"))
DimPlot(subset9,reduction = "umap",label = TRUE)
FeaturePlot(subset9, features = c("stmn1a", "si:ch211-251b21.1", "nsfl1c", "rab1ab", "hsp90ab1", "psma1" ) ,min.cutoff ="q9", pt.size=0.5)
FeaturePlot(subset16, features = c("ascl1a", "fabp7a", "ccnd1",  "si:ch211-251b21.1" ) ,min.cutoff ="q9", pt.size=0.5)

WT_3_new_subset_BP1 <- subset(WT_3_new, idents = c("BP-1"))
DimPlot(WT_3_new_subset_BP1,reduction = "umap",label = TRUE)

subset21 <-subset(WT_3_new,idents = c("RPC" ))
DimPlot(subset21,reduction = "umap",label = TRUE, pt.size=1.0)
FeaturePlot(subset21, features = c("ccnd1","stmn1a","dla", "nusap1", "foxn4", "cldn5b" ) ,min.cutoff ="q9", pt.size=1.0)




#split into positve and negative of CCND1
#pos_ids <- WhichCells(subset12, slot = 'counts', expression = ACE2 > 0 | TMPRSS2 > 0)
pos_ids <- WhichCells(subset21, slot = 'counts', expression = ccnd1 > 0)
pos_cells = subset(subset21, cells = pos_ids)
FeaturePlot(pos_cells,features = c("ccnd1","stmn1a","insm1a", "aurkb", "si:ch211-251b21.1", "fabp7a"))
FeaturePlot(pos_cells,"ccnd1")
#WhichCells(orga_SB_with_NORM, slot = 'counts', expression = ACE2 > 0)

#neg_ids <- WhichCells(Morizane_kideny_tissue, slot = 'counts', expression = ACE2 == 0 & TMPRSS2 == 0)
neg_ids <- WhichCells(subset21, slot = 'counts', expression = ccnd1 == 0)
neg_cells = subset(subset21,cells=neg_ids)
FeaturePlot(neg_cells,features = c("ccnd1","stmn1a","insm1a", "aurkb", "si:ch211-251b21.1", "fabp7a"))

pos_neg.markers <- FindMarkers(subset21, ident.1 = pos_ids, ident.2 = neg_ids, min.pct = 0.25)
head(pos_neg.markers, n = 25)
write.csv(pos_neg.markers, "subset21 WT_pos_neg_DEG.csv")

#subset method 2 hc

subset_HC <- subset(WT_3_new,idents = c("HC") )

DimPlot(subset_HC, reduction = "umap", label = FALSE, repel = TRUE, split.by = "orig.ident")

DefaultAssay(subset_HC) <- "RNA"

Subset_Cells.list <- SplitObject(subset_HC, split.by = "orig.ident")

for (i in 1:length(Subset_Cells.list)) {
  
  Subset_Cells.list[[i]] <- NormalizeData(Subset_Cells.list[[i]], verbose = FALSE)
  
  Subset_Cells.list[[i]] <- FindVariableFeatures(Subset_Cells.list[[i]], selection.method = "vst", nfeatures = 2000,
                                                 
                                                 verbose = FALSE)
  
}

reference.list <- Subset_Cells.list[c("WT_3_new3k", "p23h.new.v3")]

subset_HC.Integrated <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)

subset_HC.Integrated <- IntegrateData(anchorset = subset_HC.Integrated, dims = 1:20)

subset_HC.Integrated <- ScaleData(subset_HC, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

subset_HC.Integrated <- RunPCA(subset_HC.Integrated, npcs = 30, verbose = FALSE)

subset_HC.Integrated <- RunTSNE(subset_HC.Integrated, reduction = "pca", dims = 1:20)

subset_HC.Integrated <- FindNeighbors(subset_HC.Integrated, reduction = "pca", dims = 1:20)

subset_HC.Integrated <- FindClusters(subset_HC.Integrated, resolution = 0.5)

p1 <- DimPlot(subset_HC.Integrated, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(subset_HC.Integrated, reduction = "tsne", label = TRUE, repel = TRUE, split.by = "orig.ident") 
p2
p1
plot_grid(p1, p2)
subset_HC.markers <- FindAllMarkers(subset_HC.Integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(subset_HC.markers,file = "C:/Users/abisa/Desktop/10-WT-P23H-V3 INTEGRATED/mt-20- wt/subset_HC1.Integrated.markers.csv")
subset_HC <- FindMarkers( subset_HC.Integrated, ident.1 = "WT_3_new3k", ident.2 = "p23h.new.v3", group.by = "orig.ident"  )
write.csv(subset_HC_separate,file = "C:/Users/abisa/Desktop/10-WT-P23H-V3 INTEGRATED/mt-20- wt/subset_HC.Integrated.separate.csv")
FeaturePlot(subset_HC.Integrated, features = c("aurkb", "myca" ) ,min.cutoff ="q9", pt.size=1.0,label = TRUE, split.by = "orig.ident")
FeaturePlot(subset_HC.Integrated, features = c("aurkb", "myca" ) ,min.cutoff ="q9", pt.size=1.0,label = TRUE)
top5 <- subset_HC.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
top5
DefaultAssay(subset_HC.Integrated)  <- "integrated"
DefaultAssay(subset_HC.Integrated)  <- "RNA"
DoHeatmap(subset_HC.Integrated, features = top5$gene) + NoLegend()
DotPlot(subset_HC.Integrated, features = c( "gngt2a", "clul1","arr3a","arr3b" ) )

VlnPlot(object = subset_HC.Integrated, features = c("arr3b", "foxq2", "pde6c"), slot= "data", split.by = "orig.ident", pt.size = 0)



