# source activate r_cistopic

#library(episcanpy2r)
library(cisTopic)
library(Matrix)
library(stringr)

set.seed(2020)

start_time <- Sys.time()
print(start_time)

INPUTFILE="dgCMatrix-whole_Cusanovich_raw.rds"
dm = readRDS(INPUTFILE)

#dm = readMM(INPUTFILE)

VARFILE="ATAC_mtx_cus_whole_atlas_var.csv"
var = read.table(VARFILE, header = FALSE, sep = ",", quote = "")

OBSFILE="ATAC_mtx_cus_whole_atlas_obs_all.csv"
obs = read.table(OBSFILE, header = FALSE, sep = ",", quote = "")

#METAFILE="/home/icb/chaichoompu/Group/workspace/Kris_Workspace/scATAC-benchmarking/Real_Data/Cusanovich_2018/input/metadata.tsv"
#metadata <- read.table(METAFILE,
#                       header = TRUE,
#                       stringsAsFactors=FALSE,quote="",row.names=1)
#head(metadata)

metadata2 = as.data.frame(obs[,5])
rownames(metadata2) = obs[,2]
colnames(metadata2) = c("cell_type")
metadata = metadata2
head(metadata)

new_col = var[,1]

#reformat col_name to cisTopic format, e.g. "chr1:100-500"
if (length(grep('chr1_',new_col)) != 0){
    new_col <- str_replace(new_col,'_',':')
}
if (length(grep('chr1-',new_col)) != 0){
    new_col <- str_replace(new_col,'-',':')
}
if (length(grep('_',new_col)) != 0){
    new_col <- str_replace(new_col,'_','-')
}

colnames(dm) = obs[,2]
rownames(dm) = new_col

#dm = t(dm)

#use_condaenv("episcanpy_09102020")


#adata <- read_h5ad(INPUTFILE, conda = "episcanpy")
#dm = h5ad_to_cistopic(adata, colname = 1, rowname = 1)

cisTopicObject = createcisTopicObject(count.matrix = dm, project.name = "cusanovich_2018_seurat", keepCountsMatrix = FALSE)

#cisTopicObject <- renameCells(cisTopicObject, cellnames)

#cisTopicObject <- runCGSModels(cisTopicObject, topic=c(10, 20, 30, 40, 50, 60), seed=2020, nCores=8, burnin = 120, iterations = 150, addModels=FALSE)

######
#to test
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(10, 20, 30, 40, 50, 60, 70, 80), seed=2020, nCores=8, burnin = 120, iterations = 150, addModels=FALSE)
######

cisTopicObject <- selectModel(cisTopicObject)

logLikelihoodByIter(cisTopicObject, select=c(10, 20, 30, 40, 50, 60, 70, 80))

cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
dim(cellassign)
cellassign[1:5,1:5]

sum(colnames(cellassign) == rownames(metadata))

#saveRDS(cellassign, file = '/localscratch/chaichoompu/Cusanovich_2018/cisTopic/cisTopicObject_feature_matrices_seurat.rds')
write.csv(t(cellassign), "/localscratch/chaichoompu/Cusanovich_2018/cisTopic/cisTopicObject_feature_matrices_seurat_80topics.csv")

cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = metadata)

print("Before the first saving point")
end_time <- Sys.time()
print(end_time)
end_time - start_time

save(cisTopicObject, file = 'cisTopicObject_analysis_part_seurat_80topics.RData')

print("After the first saving point")
end_time <- Sys.time()
print(end_time)
end_time - start_time

cisTopicObject <- runUmap(cisTopicObject, target='cell')

print("Before the second saving point")
end_time <- Sys.time()
print(end_time)
end_time - start_time

save(cisTopicObject, file = 'cisTopicObject_analysis_part_seurat_80topics.RData')

print("After the second saving point")
end_time <- Sys.time()
print(end_time)
end_time - start_time

pdf("cisTopic_Cusanovich_2018_seurat_80topics.pdf")
options(repr.plot.width=10, repr.plot.height=5)
par(mfrow=c(1,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('label','pct_ReadsInPeaks'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=20)
dev.off()

print("After the last point")
end_time <- Sys.time()
print(end_time)
end_time - start_time


##from internet
## very simple export - in triplet format - to text file:
#data(CAex)
#s.CA <- summary(CAex)
#s.CA # shows  (i, j, x)  [columns of a data frame]
#message("writing to ", outf <- tempfile())
#write.table(s.CA, file = outf, row.names=FALSE)
## and read it back -- showing off  sparseMatrix():
#str(dd <- read.table(outf, header=TRUE))
## has columns (i, j, x) -> we can use via do.call() as arguments to sparseMatrix():
#mm <- do.call(sparseMatrix, dd)
#stopifnot(all.equal(mm, CAex, tolerance=1e-15))




