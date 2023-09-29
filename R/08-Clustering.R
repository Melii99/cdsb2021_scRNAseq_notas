########## Clustering ##########

### Dataset ilustrativo: 10X PBMC4k no filtrado ###

## Dataset de células mononucleares de sangre periférica humana (PBMC) de 10X Genomics


## Descarga de datos  ##
library(BiocFileCache)

bfc <- BiocFileCache()

raw.path <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples",
  "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
))

untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
library(Matrix)

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)


## Anotación ##

library(scater)

rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol
)

library(EnsDb.Hsapiens.v86)

## Mapear genes
location <- mapIds(EnsDb.Hsapiens.v86,
                   keys = rowData(sce.pbmc)$ID,
                   column = "SEQNAME", keytype = "GENEID"
)


##  Detección de empty droplets  ##

## Filtro de células, droplets vacios (FDR .001)
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]


## Control de calidad QC ##

## QC metrics de genes  mitocondriales
stats <- perCellQCMetrics(sce.pbmc,
                          subsets = list(Mito = which(location == "MT"))
)
## Obtener genes mitocondriales con alta expresión
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type = "higher"
)
## Filtrar genes mitocondriales con alta expresión
sce.pbmc <- sce.pbmc[, !high.mito]


## Normalización ##

## Normalización por quick clusters
library(scran)

set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)

## Transformación logarítmica
sce.pbmc <- logNormCounts(sce.pbmc)


## Modelado de la varianza ##

## Modelado de la varianza usando como supuesto una distribución de Poisson (no ERCC)
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)

## Obtenemos HVGs
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1) # Los top genes no se guardaron en sce,pbmc


## Reducción de dimensionalidad ##

set.seed(10000)

## PCA usando los HVGs
sce.pbmc <- denoisePCA(sce.pbmc,
                       subset.row = top.pbmc,
                       technical = dec.pbmc # Ya se está usando "filtro" para los PCs
)

## t-SNE
set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred = "PCA")

## UMAP
set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred = "PCA")
