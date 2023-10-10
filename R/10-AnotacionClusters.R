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



### Identificando los genes con anotación dirigida ###


## ¿Por qué las células en este clúster se etiquetan como el tipo celular X?

## Examina la expresión de los genes marcadores para cada etiqueta en el dataset de prueba

## Si una célula en el dataset de prueba está asignado con confianza a una
## etiqueta en particular, uno esperaría que tenga una fuerte expresión de los
## marcadores de esa etiqueta (al menos sobreexpresión con respecto a las
## células asignadas a otras etiquetas)

# install gmp, ClusterR, mbkmeans dependencies if needed

## Agregar las etiquetas de las células predichas por SingleR a sce.pbmc
sce.pbmc$labels <- pred$labels

## Obtener los genes marcadores de las metadatos (de.genes) del objeto pred almacenar en all.markers.
all.markers <- metadata(pred)$de.genes

## Ejemplo de label: células B
lab <- "B-cells"


## Obtener los 10 genes principales de los genes marcadores específicos para las
## células tipo B y almacenarlos en top.markers.
top.markers <- Reduce(union, sapply(all.markers[[lab]], head, 10)) # head, 10

## Heatmap de la expresión de los genes marcadores específicos (top.markers) para las células tipo B
plotHeatmap(sce.pbmc,
            order_columns_by = "labels",
            features = top.markers, center = TRUE, zlim = c(-3, 3), main = lab
)



### Comparando las etiquetas con los clústeres ###


## Queremos saber cuántas células de cada clúster han sido asignadas a cada etiqueta
## predicha por SingleR (anotaciones de las referencias)

##  Generar una tabla que muestra la relación entre las etiquetas asignadas por
## SingleR y los clústeres que encontramos y guardamos en  sce.pbmc
tab <- table(Assigned = pred$pruned.labels, Cluster = sce.pbmc$cluster)


## Heatmap de la proporción de células en cada clúster que han sido asignadas a cada etiqueta
library(pheatmap)
pheatmap(prop.table(tab, margin = 2),
         color = colorRampPalette(c("white", "blue"))(101)
)

## Heatmap simlar al anterior, ahora usando el logaritmo base 2 del número de
## células en cada clúster asignadas a cada etiqueta predicha por SingleR
## Se añade un pseudo-recuento de 10 para evitar saltos de color abruptos con solo 1 célula
pheatmap(log2(tab + 10),
         color = colorRampPalette(c("white", "blue"))(101)
)


## Visualizando ##

## t-SNE donde las células se encuentran coloreadas según las etiquetas dadas por SingleR
plotTSNE(sce.pbmc, colour_by = "labels", text_by = "labels")

##  t-SNE donde las células se encuentran coloreadas según los clústeres encontrados
plotTSNE(sce.pbmc, colour_by = "cluster", text_by = "labels")


### Resumen de la anotación basada en una referencia (e.g., SingleR) ###

##  centra en aspectos de los datos que se sabe son interesantes, simplifica el
## proceso de la interpretación biológica

## Está restringido por la diversidad y la resolución de las etiquetas disponibles
## en el dataset de referencia

## También, es posible suplir referencias personalizadas a SingleR



### Asignando las etiquetas de tipos celulares a partir de marcadores ###

## ¿Para qué usar nuestros genes marcadores agrupados?

# Revisarlos en hojas de cálculo
# Observar heatmaps
# Realizar un gene set enrichment analysis


## Gene set enrichment analysis ##

## Identifica las rutas y procesos que están (relativamente) activos en cada clúster
## basado en la sobreexpresión de los genes asociados en comparación con otros clústeres

## Método confiable para determinar si las rutas están sobre- o sub- expresadas entre clúesteres

## Existen un montón de herramientas para gene set enrichment analysis

## Todas las conclusiones son relativas a otros clústeres, haciéndo más difícil
## determinar la identidad celular si alguno no está presente en el mismo estudio



## Calculando las actividades de los conjuntos de genes ###

## Calcular el promedio de la expresión en log en todos los genes, en un conjunto
## de genes para cada célula y examinar los clústeres con valores altos (gene set activities)

## Se necesita proveer de conjuntos de genes

##  todos los genes en el conjunto pueden exhibir el mismo patrón de diferencia
## y los genes no-DE añadirán ruido, “diluyendo” la fuerza de cualquiera de las
## diferencias comparadas a un análisis que se centra directamente en genes DE

## Es más una visualización útil que la base para cualquier análisis estadístico real



### Resumen y recomendaciones ###

## La anotación de tipos celulares “automática”, como SingleR, es mejor cuando funciona
## (i.e. cuando hay un dataset de referencia apropiado)

## Usualmente necesitaremos usar un método manual, como aquellos basados en agrupar
## los genes marcadores (e.g., gene set enrichment analysis)

## La anotación del tipo celular ofrecerá una reconsideración inmediata de los
## parámetros del agrupamiento y/o algunos retoques manuales a los clústeres



### Información de la sesión de R ###
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
