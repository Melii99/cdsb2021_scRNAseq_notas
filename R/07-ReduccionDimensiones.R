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

## Modelado de la varianza ##

set.seed(1001)
## Modelado de la b¿varianza usando el supuesto de una distribución poisson
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)

## Obtenemos el  10% de los genes con mayor variación en expresión como HVGs
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1)



### Análisis de Componentes Principales ###


## Es el "arma principal" de la reducción de dimensionalidad
## Descubre las combinaciones (lineales) de “features” que capturan la cantidad más grande de variación

## En un PCA, la primer combinación lineal (componente principal) se elige tal que
## permite capturar la mayor varianza a través de las células.
## El siguiente PC se elige tal que es “ortogonal” al primero y captura la cantidad
## más grande de la variación restante, y así sucesivamente


## PCA aplicado a datos de scRNA-seq ##

## Podemos realizar reducción de dimensionalidad al aplicar PCA en la matriz de cuentas transformadas
## (log-counts matrix) y restringiendo los análisis posteriores a los primeros PCs (top PCs)

## PROS:
## Reducción de dimensiones sin perder demasiada información
## La técnica de PCA tiene muchas propiedades teóricas bien estudiadas.
## Hay varias formas rápidas de realizar PCA en datasets grandes.


## Suposiciones de PCA aplicadas a los datos de scRNA-seq ##

## Los procesos biológicos afectan múltiples genes en una manera coordinada

## Los primeros PCs probablemente representan la estructura biológica dado que más
## variación puede ser capturada considerando el comportamiento correlacionado de muchos genes

## Se espera que el ruido técnico azaroso afecte cada gen independientemente

## Consideración: Los primeros PCs capturarón “batch effects” (efectos de lote)
## que afectan muchos genes en una manera coordinada


## PCA con runPCA ##

## Por default, runPCA() usa un método rápido aproximado que realiza simulaciones,
##por lo tanto, es necesario ‘configurar la semilla’ para obtener resultados reproducibles

## PCA de los topGenes obtenidos para zeisel
library(scran)
library(scater)

set.seed(100)
sce.zeisel <- runPCA(sce.zeisel,
                     subset_row = top.zeisel
)



### Eligiendo el número de PCs ###


## Esta elección en análoga a la elección del numero de HVGs.
## Elegir más PCs evitará descartar señal biológica a expensas de retener más ruido

## Es común seleccionar un número de PCs “razonable” pero arbitrario (10-50),
## continuar con el análisis y regresar para checar la robustez de los resultados


## Estrategias data-driven para seleccionar el número de PCs ##

## Método del codo ##

## Una heurística simple es elegir el número de PCs basado en el porcentaje de
## varianza explicado por PCs sucesivos

library(PCAtools)

percent.var <- attr(reducedDim(sce.zeisel), "percentVar") # reducedDim()
chosen.elbow <- PCAtools::findElbowPoint(percent.var) # findElbowPoint()

## Graficar el porcentaje de la varianza explicada
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")


## Basados en la estructura de la población ##

## Esta es una aproximación heurística más sofisticada que usa el número de
## cluster como un proxy del número de subpoblaciones

## Supongamos que esperamos d subpoblaciones de células, en ese caso, necesitamos
## d-1 dimensiones para garantizar la separación de todas las subpoblaciones

## Pero en un escenario real realmente no sabes cuántas poblaciones hay !!!

## Intentae con un rango para d y únicamente considera valores que produzcan a lo más d+1 clusters
## Cuando se seleccionan más clusters con menos dimensiones se produce ‘overclustering’
## Elegir una d que maximice el número de clusters sin caer en ‘overclustering’


## PROS:
## Es una solución pragmática que soluciona el equilibrio sesgo-varianza en los
## análisis posteriores (especialmente clustering)

## CONTRAS:
## Hace suposiciones fuertes sobre la naturaleza de las diferencias biológicas entre
##los clusters, y de hecho supone la existencia de clusters, los cuales podrían no
## existir en procesos biológicos como la diferenciación


choices <- getClusteredPCs(reducedDim(sce.zeisel))
chosen.clusters <- metadata(choices)$chosen

plot(choices$n.pcs, choices$n.clusters,
     xlab = "Number of PCs", ylab = "Number of clusters"
)
abline(a = 1, b = 1, col = "red")
abline(v = chosen.clusters, col = "grey80", lty = 2)
