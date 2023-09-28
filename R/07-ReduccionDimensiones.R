########## Reducción de dimensiones ##########

## El siguiente paso en el análisis de scRNA-seq usualmente consiste en identificar
## grupos de células “similares”

## Por ejemplo: un análisis de clustering busca identificar células con un perfil
## transcriptómico similar al calcular distancias entre ellas

## Si tuviéramos un dataset con dos genes podríamos hacer una gráfica de dos dimensiones
## para identificar clusters de células. Pero… tenemos decenas de miles de genes.
## Por ello es relevante la reducción de dimensionalidad

## La expresión de diferentes genes estaría correlacionada si estos genes son
## afectados por el mismo proceso biológico !!!

## PROS:
## Reduce el trabajo computacional en análisis posteriores
## Reduce el ruido al “promediar” mútiples genes obteniendo una representación
## más precisa de los patrones en los datos
## Permite una graficación efectiva en dos dimensiones



### Dataset ilustrativo: Zeisel ###

## Estudio de tipos celulares en el cerebro de ratón (oligodendrocitos, microglia,
## neuronas, etc) procesados con el sistema STRT-seq (similar a CEL-Seq)

## DEscarga de los datos
library(scRNAseq)
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)

## Estos datos contienen tipos celulares previamente anotados
table(sce.zeisel$level1class)

## Control de calidad QC ##

## Descartar celulas con alto contenido mitocondrial o con alto porcentaje de spike-ins
library(scater)

is.mito <- which(rowData(sce.zeisel)$featureType == "mito")

stats <- perCellQCMetrics(sce.zeisel,
                          subsets = list(Mt = is.mito)
)

qc <- quickPerCellQC(stats,
                     percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent")
)

sce.zeisel <- sce.zeisel[, !qc$discard]

## Normalizacion ##

## Encontrar  clusters rápidos para las células y usar esa información para calcular los factores de tamaño
library(scran)

set.seed(1000)
clusters <- quickCluster(sce.zeisel)

sce.zeisel <- computeSumFactors(sce.zeisel,
                                cluster = clusters
)

sce.zeisel <- logNormCounts(sce.zeisel)


## Modelado de la varianza ##

## Modelar la varianza a partir de los ERCC
dec.zeisel <- modelGeneVarWithSpikes(sce.zeisel, "ERCC")

## Obtener los 2000 genes "más relen¿vantes" con mayor varianza biológica
top.zeisel <- getTopHVGs(dec.zeisel, n = 2000)



### Dataset ilustrativo: 10x PBMC4k no filtradas ###

## Dataset “Células mononucleares humanas de sangre periférica” de 10X Genomics

## Descarga de los datos
library(BiocFileCache)

bfc <- BiocFileCache()

raw.path <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples",
  "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
))

untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

## Leer en R
library(DropletUtils)
library(Matrix)

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

## Anotación ##

## Anotación de los genes
library(scater)

rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol
)

library(EnsDb.Hsapiens.v86)

location <- mapIds(EnsDb.Hsapiens.v86,
                   keys = rowData(sce.pbmc)$ID,
                   column = "SEQNAME", keytype = "GENEID"
)

## Detección de droplets vacíos ##

# Detección de droplets vacíos con enptyDrops
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))

## Selección de droplets con células (con un false discovery rate igual o menor a .001)
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]


## Control de calidad QC ##

## Identificación de genes mitocondriales
stats <- perCellQCMetrics(sce.pbmc,
                          subsets = list(Mito = which(location == "MT"))
)

## Genes mitocondriales altamente expresados
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type = "higher"
)

## Eliminar genes mitocondriales altamente expresados
sce.pbmc <- sce.pbmc[, !high.mito]


## Normalización ##

## Normalización por quick clusters y obtención de sizefactors
library(scran)

set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)

## Transformación logarítmica
sce.pbmc <- logNormCounts(sce.pbmc)


