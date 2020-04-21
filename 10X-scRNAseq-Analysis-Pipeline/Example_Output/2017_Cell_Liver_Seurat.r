sink("2017_Cell_Liver_Seurat/logofcode.txt")
library(Seurat)
library(dplyr)
library(Matrix)
single.data0 <- Read10X(data.dir = "CellRangerFolder/")
at_least_one0 <- apply(single.data0, 2, function(x) sum(x>0))
pdf("2017_Cell_Liver_Seurat/1_2017_Cell_Liver_single.data0_Histogram_At_Least_One.pdf")
hist(at_least_one0, breaks = 100,
main = "Distribution of detected genes",
xlab = "Genes with at least one tag")
dev.off()

pdf("2017_Cell_Liver_Seurat/2_2017_Cell_Liver_Single.Data0_Histogram_Expr_Sum_per_Cell.pdf")
hist(colSums(single.data0),
breaks = 100, main = "Expression sum per cell",
xlab = "Sum expression")
dev.off()

tmp <- apply(single.data0, 1, function(x) sum(x>0))
table(tmp>=3)

keep <- tmp>=3
tmp <- single.data0[keep,]
at_least_one0 <- apply(tmp, 2, function(x) sum(x>0))
summary(at_least_one0)

dim(tmp)

single.data0 <- CreateSeuratObject(counts = single.data0, min.cells = 3, min.features = 200)
pdf("2017_Cell_Liver_Seurat/4_2017_Cell_Liver_data0Total_Expr_Before_Norm.pdf")
hist(colSums(single.data0$RNA@data),
breaks = 100,
main = "Total expression before normalisation",
xlab = "Sum of expression")
dev.off();
single.data0 <- NormalizeData(object = single.data0, normalization.method = "LogNormalize", scale.factor = 1e4)
pdf("2017_Cell_Liver_Seurat/5_2017_Cell_Liver_data0_Total_Expr_After_Norm.pdf")
hist(colSums(single.data0$RNA@data),
breaks = 100,
main = "Total expression after normalisation",
xlab = "Sum of expression")
dev.off();
single.data0 <- FindVariableFeatures(object = single.data0, selection.method = "vst")
single = single.data0
pdf("2017_Cell_Liver_Seurat/4_2017_Cell_Liver_Total_Expr_Before_Norm.pdf")
hist(colSums(single$RNA@data),
breaks = 100,
main = "Total expression before normalisation",
xlab = "Sum of expression")
dev.off();
single <- NormalizeData(object = single, normalization.method = "LogNormalize", scale.factor = 1e4)
pdf("2017_Cell_Liver_Seurat/5_2017_Cell_Liver_Total_Expr_After_Norm.pdf")
hist(colSums(single$RNA@data),
breaks = 100,
main = "Total expression after normalisation",
xlab = "Sum of expression")
dev.off();
single <- FindVariableFeatures(object = single, mean.function = "FastExpMean", dispersion.function = "FastLogVMR", mean.cutoff=c(0.0125, 3), dispersion.cutoff = c(0, 0.5))
single <- ScaleData(object = single)
single <- RunPCA(object = single, pc.genes = single@var.genes)
single <- FindNeighbors(object = single)
single <- FindClusters(object = single, reduction.type = "pca",dims.use = 1:10, resolution = 0.6)
single <- RunTSNE(object = single, check_duplicates = FALSE, dims.use = 1:10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.6)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.6_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.6_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.6_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.6)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.6_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.6_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.6_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.6)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.6_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.6_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.6_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.6)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.6_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.6_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.6_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.6)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.6_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.6_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.6_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.6)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.6_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.6_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.6_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.4)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.4_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.4_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.4_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.4)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.4_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.4_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.4_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.4)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.4_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.4_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.4_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.4)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.4_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.4_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.4_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.4)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.4_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.4_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.4_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.4)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.4_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.4_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.4_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.9)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.9_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.9_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.9_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.9)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.9_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.9_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.9_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.9)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.9_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.9_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.9_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.9)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.9_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.9_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.9_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.9)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.9_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.9_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.9_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.9)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.9_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.9_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.9_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.1)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.1_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.1_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.1_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.1)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.1_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.1_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.1_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.1)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.1_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.1_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.1_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.1)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.1_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.1_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.1_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.1)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.1_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.1_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.1_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.1)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.1_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.1_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.1_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.2)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.2_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.2_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.2_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.2)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.2_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.2_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.2_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.2)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.2_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.2_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.2_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.2)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.2_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.2_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.2_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.2)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.2_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.2_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.2_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.2)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.2_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.2_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.2_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.3)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.3_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.3_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.3_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.3)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.3_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.3_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.3_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.3)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.3_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.3_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.3_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.3)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.3_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.3_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.3_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.3)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.3_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.3_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.3_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.3)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.3_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.3_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.3_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.5)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.5_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.5_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.5_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.5)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.5_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.5_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.5_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.5)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.5_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.5_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.5_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.5)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.5_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.5_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.5_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.5)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.5_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.5_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.5_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.5)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.5_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.5_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.5_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.7)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.7_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.7_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.7_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.7)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.7_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.7_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.7_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.7)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.7_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.7_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.7_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.7)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.7_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.7_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.7_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.7)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.7_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.7_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.7_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.7)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.7_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.7_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.7_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.8)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.8_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.8_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.8_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.8)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.8_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.8_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.8_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.8)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.8_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.8_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.8_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.8)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.8_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.8_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.8_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.8)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.8_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.8_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.8_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 0.8)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res0.8_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res0.8_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res0.8_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 1.0)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 30)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res1.0_Per30.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res1.0_Per30.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res1.0_Per30.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 1.0)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 50)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res1.0_Per50.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res1.0_Per50.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res1.0_Per50.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 1.0)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 100)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res1.0_Per100.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res1.0_Per100.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res1.0_Per100.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 1.0)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 10)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res1.0_Per10.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res1.0_Per10.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res1.0_Per10.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 1.0)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 20)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res1.0_Per20.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res1.0_Per20.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res1.0_Per20.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();
single <- FindClusters(object = single, dims.use = 1:10, resolution = 1.0)
single <- RunTSNE(object = single, dims.use = 1:10, check_duplicates = FALSE, perplexity = 40)
write.table(single@active.ident, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_CellsIdentity_Res1.0_Per40.txt",sep="\t")
write.table(single$tsne@cell.embeddings, file = "2017_Cell_Liver_Seurat/2017_Cell_Liver_TSNE_Res1.0_Per40.txt",sep="\t")

pdf("2017_Cell_Liver_Seurat/6_2017_Cell_Liver_PCAPlot.pdf")
DimPlot(object = single, reduction = "pca")
dev.off();

pdf("2017_Cell_Liver_Seurat/8_2017_Cell_Liver_TSNEPlot_Res1.0_Per40.pdf")
DimPlot(object = single, reduction = "tsne")
dev.off();

