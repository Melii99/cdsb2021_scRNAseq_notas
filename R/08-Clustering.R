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



### Clustering ###

## Es un procedimiento no supervisado par definir grupos de células con perfiles de expresión similares
## Su propósito principal es resumir los datos en un formato digerido susceptible a interpretación humana
## Nos permite asignar etiquetas (por ejemplo, tipos celulares) a las células

## ¿Por qué no realizamos el clustering sobre las coordenadas de t-SNE/UMAP?

## Las técnicas de t-SNE/UMAP han comprimido datos altamente multi-dimensionales en dos dimensiones
## Esta compresión inevitablemente ha provocado la perdida de información
## Por lo tanto, agrupamos sobre los PCs y después visualizamos las identidades de
## los clusters en la gráfica t-SNE/UMAP


## ¿Cuál es el verdadero clustering? ##

## Un cluster no implica un tipo celular !!!
## Podemos definir tantos clusters como queramos y podemos utilizar el algoritmo que más nos acomode
## El clustering, como un microscopio, simplemente es una herramienta para explorar los datos



### Clustering basado en grafos ###

## El clustering basado en grafos fue popularizado (más no inventado) por su uso en Seurat
## Tiene como objetivo construir un grafo en el que cada nodo es una célula que está
## conectada a sus vecinos más cercanos en el espacio multidimensional

## Vecinos más cercanos; aquellas otras células que tienen perfiles de expresión muy similares


## Gráficas KNN - k-nearest neighbors ##

## Preguntamos para cada célula cuáles son los k de vecinos más cercanos
## (Células con expresión más semejante)

## Ej. si tengo 20 células y busco k = 19 (20 - 1) vamos a encontrar a todas las células conectadas

## El parámetro 'k' es lo que ayuda a resolver la granularidad o resolución de los clusters
## (qué tan agrupados y grandes los queremos, a menor k, más pequeños y más resolución)


## Gráficas SNN - Shared nearest neighbors ##

## A partir de una gráfica KNN se puede construir una grafica SNN

## En este tipo de grafo, dos células estarán conectadas por una arista si son vecinos
## más próximos o si comparten alguno de sus vecinos más próximos.

## También podemos asignar pesos a cada arista del grafo, basándonos en la similaridad de las células
## involucradas, dándole pesos más altos a células que están más cercanamente relacionadas


## Gráficas SNN con pesos en las aristas ##

## Cómo asignar los pesos a las arístas:

## Rango: El peso entre dos nodos está dado por k-r/2 donde r es la suma más pequeña
## de los rangos (de proximidad, el vecino más cercano tiene el rango 1) para cualquiera
## de los vecinos compartidos (vecino más cercano = 1, el que sigue en proximidad = 2...)

## Número: el peso entre dos nodos es igual al número de vecinos más próximos compartidos

## Jaccard: el peso entre dos nodos es igual a la similaridad de Jaccard entre los
## conjuntos de vecinos de estos nodos ()


## Obteniendo comunidades a partir de una gráfica SNN pesada mediante un algoritmo de clustering ##

## A partir de una gráfica SNN pesada podemos aplicar algoritmos para identificar comunidades de células

## Comunidades son grupos de células que están más conectadas a células en el mismo
## grupo que lo que están a células de un grupo diferente

## Cada comunidad representa un cluster



### Resumen de clustering basado en grafos ###

## La construcción y búsqueda de una red KNN es rápida, por lo tanto, es escalable para datasets grandes

## Debes evitar obtener conclusiones fuertes acerca de la forma de los clusters o
## la distribución de células dentro de cada cluster

## El algoritmo, conecta cada célula con un número mínimo de células vecinas, lo cual
## reduce el riesgo de clusters no informativos con unos pocos outliers

## Después de la construcción del grafo, no se almacena información adicional más alla de las
## células vecinas. Esto puede producir subclusters artificiales en regiones con muchas células !!!

## CONSIDERACIONES:
## ¿Cuántas céulas vecinas debo considerar durante la construcción del grafo?
## ¿Cómo debo pesar las aristas?
## ¿Cuál algoritmo de detección de comunidades se debe usar para definir los clusters?



### Implementación de clustering basado en grafos ###


## Construcción de un grafo usando una k = 10 (usar a los 10 vecinos mpas cercanos)
library(scran)
g <- buildSNNGraph(sce.pbmc, k = 10, use.dimred = "PCA")

## Identificar las comunidades (clusters) usando el método Walktrap
clust <- igraph::cluster_walktrap(g)$membership

## Visualización de los clusters con un plot t-SNE
library(scater)
sce.pbmc$cluster <- factor(clust)
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster")


## ¿Qué pasa si utilizas una k más grande o más pequeña? ##


## Construcción de un grafo usando una k = 50 (usar a los 50 vecinos mpas cercanos)
library(scran)
g50 <- buildSNNGraph(sce.pbmc, k = 50, use.dimred = "PCA")

## Identificar las comunidades (clusters) usando el método Walktrap
clust50 <- igraph::cluster_walktrap(g50)$membership

## Visualización de los clusters con un plot t-SNE
library(scater)
sce.pbmc$cluster50 <- factor(clust50)
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster50")


## En esta implementación:
## La construcción de la red KNN se baso en la distancia Euclideana entre células
## La construcción de la red KNN implica que las aristas se crean entre todos los
## pares de células que comparten por lo menos un vecino
