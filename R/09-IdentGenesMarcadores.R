########## Identificación de Genes Marcadores ##########


## Ahora que hemos obtenido los clústeres, nos preguntamos, pero qué son?
## (e.g. ¿qué tipo celular es el clúster 1?)

## ¿Cuáles genes están dirigiendo el agrupamiento (e.g., ¿cuáles son los genes
## diferencialmente expresados entre los clústeres 1 y 2?)

## Idea: Mirar las diferencias en los perfiles de expresión de las células de los diferentes clústeres


### Dataset ilustrativo: PBMC4k 10X sin filtrar ###

## Dataset “Células mononucleares humanas de sangre periférica” de 10X Genomics

## Descarga de datos ##

# Descarga de datos de  pbmc4k
library(BiocFileCache)

bfc <- BiocFileCache()

raw.path <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples",
  "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
))

untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

## Leer datos en R
library(DropletUtils)
library(Matrix)

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)


## Anotación ##

# Asignación de 'etiquetas' de los genes
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol
)

## Mapear los genes
library(EnsDb.Hsapiens.v86)

location <- mapIds(EnsDb.Hsapiens.v86,
                   keys = rowData(sce.pbmc)$ID,
                   column = "SEQNAME", keytype = "GENEID"
)

## Control de calidad QC ##

## Detección de droplets vacíos y filtrado (por false discovery rate menor a 0.001)
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

## Identificar y obtener métricas QC de genes mitocondriales
stats <- perCellQCMetrics(sce.pbmc,
                          subsets = list(Mito = which(location == "MT"))
)

## Obtener genes mitocondriales con alta expresión
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type = "higher"
)
## Filtrar los genes mitocondriales con alta expresión
sce.pbmc <- sce.pbmc[, !high.mito]


## Normalización ##

library(scran)
set.seed(1000)

## Normalización por quickclusters
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)

## Transformación logaritmica de la matríz de cuentas
sce.pbmc <- logNormCounts(sce.pbmc)


## Genes variables HGVs ##

## Identificación de los genes altamente variables HVGs usando como supuesto una distrib. poisson
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1)


## Reducción de dimensiones ##

## Reducción de dimensiones por PCA
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc,
                       subset.row = top.pbmc, # Solo se utilizarán llos HVGs
                       technical = dec.pbmc # Especifica la información de variabilidad genética modelada
)

## tSNE para visualización
set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred = "PCA")
## UMAP para visualización
set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred = "PCA")

## Clustering ##

## Construcción de un grafo SNN
g <- buildSNNGraph(sce.pbmc, k = 10, use.dimred = "PCA")

## Identificación de los clusters a partir del grafo usando el algoritmo walktrap
clust <- igraph::cluster_walktrap(g)$membership

## Agregar la información de los clusters en el sce.pbmc
sce.pbmc$cluster <- factor(clust)
