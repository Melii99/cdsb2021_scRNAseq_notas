########## Anotación de clusters de células ##########


## Ahora estamos a punto de obtener la interpretación biológica de los resultados
## Esta es la tarea más retadora en los análisis de datos scRNA-seq

## La obtención de clústeres es más o menos directa
## ¿Cuál es el estado biológico que está representado por cada uno de los clústeres?

## Necesitamos hacer un puente entre el gap del dataset actual y el conocimiento
## biológico a priori (no siempre está disponible en una forma consistente y cualitativa)

## ¿Qué es un tipo celular?


## Aplicaremos varios métodos computacionales que explotan la información a priori
## para asignar el significado a un dataset no caracterizado de scRNA-seq.


## Algunas fuentes de información a priori ##

## Conjuntos de genes curados (e.g. Gene Ontology)
## Perfiles de expresión de bases de datos publicadas de referencia
## Los datos raros que tú hayas escondido en tu cerebro
## Google



### Dataset ilustrativo: PBMC4k 10X sin filtrar ###

## Dataset “Células mononucleares humanas de sangre periférica” de 10X Genomics

## Descarga de datos ##

## Descarga de datos de pbmc4k
library(BiocFileCache)

bfc <- BiocFileCache()

raw.path <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples",
  "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
))

untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

## Leer los datos en R
library(DropletUtils)
library(Matrix)

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)


## Anotación ##

## Asignación de 'etiquetas' de los genes
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

## Normalización por quickclusters
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)

## Transformación logaritmica de la matríz de cuentas
sce.pbmc <- logNormCounts(sce.pbmc)


## Genes variables HVGs ##

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



### Continuación: Asignando las etiquetas celulares a partir de los datos de referencia ###


## Un enfoque directo es comparar los perfiles de expresión single-cell con datasets previamente anotados

## Las etiquetas pueden entonces ser asignadas a cada célula en nuestro dataset
## no caracterizado de prueba basado en la muestra de referencia más similar,
## por dar alguna definición de “similar”

## Cualquier dataset de expresión génica etiquetado (microarreglos, RNA-seq bulk,
## scRNA-seq) puede ser usado como una referencia

## Sin embargo, su confiabilidad depende enormemente en la calidad de los datos
## originales y la experiencia de los autores originales quienes asignaron las
## etiquetas en primer lugar

## Asignar las etiquetas a un dataset de “prueba” a partir de un dataset de “entrenamiento”
## (referencia), es un problema estándar en estadística / machine learning

## Es preferible tener una idea de lo que podríamos encontrar parra elegir un
## dataset de referencia

## Hay que revisar la calidad de los datos y la forma en que se procesaron

## Usaremos el método SingleR (Aran et al. 2019)



## SingleR ##

## Asigna las etiquetas a las células basado en las muestras de referencia con las
## correlaciones de rangos más altas de Spearman

## Para reducir el ruido, identifica genes marcadores entre pares de etiquetas
## (en la referencia) y calcula la correlación usando solamente esos marcadores

## Hace algún tipo de tuneado fino, repitiendo las correlaciones solamente con los
## genes marcadores de las etiquetas con el mejor score, ayudando a resolver cualquier
## ambigüedad entre esas etiquetas al eliminar el ruido a partir de marcadores
## irrelevantes para otras etiquetas


##  Referencias ##

# Human #
#celldex::BlueprintEncodeData()
#celldex::DatabaseImmuneCellExpressionData()
#celldex::HumanPrimaryCellAtlasData()
#celldex::MonacoImmuneData()
#celldex::NovershternHematopoieticData()

# Mice #
#celldex::ImmGenData()
#celldex::MouseRNASeqData()



## Usando las referencias ##

## Definimos las referencias a usar
library(celldex)
ref <- celldex::BlueprintEncodeData()


## Usando las referencias integradas ##

## Especificar las etiquetas y clasificar según las asignadas por SingleR().
library(SingleR)

pred <- SingleR(
  test = sce.pbmc, ref = ref,
  labels = ref$label.main
)

## Visualizar la asignación de etiquetas obtenidas por SingleR
## Heatmap de los scores por célula y por etiqueta
plotScoreHeatmap(pred)


## Idealmente, cada célula debería exhibir un score alto en una etiqueta
## relativa a todas las otras

## Los scores se muestran antes de cualquier tuneado fino y son normalizadas
## a [0, 1] dentro de cada célula



### Podado de etiquetas (Label pruning) ###


## SingleR intentará podar aquellas asignaciones de baja calidad marcándolas como NA

## El podado se hace calculando la diferencia del score de la etiqueta asignada
## a partir del score de la mediana dentro de cada célula y entonces podando las
## células con un valor pequeño de esta diferencia

## Se guardan las etiquetas "podadas"
total_pruned <- sum(is.na(pred$pruned.labels))

## Heatmap de los scores por célula y por etiqueta de las etiquetas "podadas"
plotScoreHeatmap(pred, show.pruned = TRUE)

## Plot de la distribución de los scores de las clasificaciones realizadas por SingleR
plotScoreDistribution(pred)