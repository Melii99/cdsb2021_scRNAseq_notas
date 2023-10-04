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


## Agregar los clusters a colData de nuestro objeto
library(scater)
sce.pbmc$cluster <- factor(clust)

## Visualización de los clusters coloreados con un plot t-SNE
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster")


## ¿Qué pasa si utilizas una k más grande o más pequeña? ##


## Construcción de un grafo usando una k = 50 (usar a los 50 vecinos mpas cercanos)
library(scran)
g50 <- buildSNNGraph(sce.pbmc, k = 50, use.dimred = "PCA")

## Identificar las comunidades (clusters) usando el método Walktrap
clust50 <- igraph::cluster_walktrap(g50)$membership

## Agregar los clusters a colData de nuestro objeto
library(scater)
sce.pbmc$cluster50 <- factor(clust50)

## Visualización de los clusters coloreados con un plot t-SNE
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster50")


## En esta implementación:
## La construcción de la red KNN se baso en la distancia Euclideana entre células
## La construcción de la red KNN implica que las aristas se crean entre todos los
## pares de células que comparten por lo menos un vecino



### Eligiendo un valor de k ###

## El valor de k puede ser toscamente interpretado como el tamaño anticipado de la
## subpoblación más pequeña

## Si una subpoblación tiene menos que (k+1) células entonces el método será forzado
## a construir aristas entre células de esa subpoblación y células de otras subpoblaciones

## Esto incrementa el riesgo de que la subpoblación en cuestión no forme su propio cluster

## MIENTRAS MÁS PEQUEÑA SEA LA K, MÁS CLUSTERS SE VAN A GENERAR !!!



### Una implementación diferente: estilo Seurat ###

# Jaccard-based weights followed by Louvain clustering
# aka 'Seurat-style' clustering

## Construcción de un grafo usando una k = 10 y calcula la similitud de Jaccard
## entre los vecinos más cercanos
g2 <- buildSNNGraph(sce.pbmc, k = 10, use.dimred = "PCA", type = "jaccard")

## Identificar las comunidades (clusters) usando el algoritmo de Louvain
clust2 <- igraph::cluster_louvain(g2)$membership

## Agregar los clusters a colData de nuestro objeto
sce.pbmc$cluster2 <- factor(clust2)

## Visualización de los clusters coloreados con un plot t-SNE
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster2")



### Detalles de las implementaciones más comunes ###

## Pipelines basados en Seurat:

# Pesos basados en Jacard
# Clustering Louvain

## Pipelines basados en Scran:

# Pesos basados en Randos
# Clustering Walktrap

## Comparando estas 2 implementaciones ##

library("patchwork")

## Grafica de las tSNEs comparado las comunidades generadas por clustering
## basado en Scran y clustering basado en Seurat
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster") +
  plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster2")


### Otras implementaciones ###

## Distintas métricas de distancia

## construye un grafo de vecinos más cercanos (SNN: Shared Nearest Neighbor)
g.num <- buildSNNGraph(sce.pbmc, use.dimred = "PCA", type = "number")

## construye un grafo SNN donde las conexiones en el grafo se basan
##  en la similitud de Jaccard entre los vecinos cercanos de las células.
g.jaccard <- buildSNNGraph(sce.pbmc, use.dimred = "PCA", type = "jaccard")

## construye un grafo de vecinos más cercanos (KNN: K-Nearest Neighbors)
g.none <- buildKNNGraph(sce.pbmc, use.dimred = "PCA")

## Distintos métodos de clustering

## Algoritmo de Louvain.Encuentra comunidades maximizando la modularidad en el grafo.
## A menudo utilizado para detectar comunidades en grafos grandes.
clust.louvain <- igraph::cluster_louvain(g)$membership

## Algoritmo Infomap. Se basa en la idea de que las comunidades se corresponden
## con las rutas óptimas para comprimir información en un flujo de información.
## A menudo utilizado para encontrar estructuras modulares en redes complejas.
clust.infomap <- igraph::cluster_infomap(g)$membership

##A lgoritmo de Greedy Fast. Iterativamente agrega nodos a comunidades para maximizar
## la modularidad. Es un algoritmo eficiente para encontrar comunidades en grafos grandes.
clust.fast <- igraph::cluster_fast_greedy(g)$membership

## Algoritmo Label Propagation. Los nodos se asignan a la comunidad de su vecino con mayor
## frecuencia. Es un algoritmo simple y rápido para detectar comunidades en grandes redes.
clust.labprop <- igraph::cluster_label_prop(g)$membership

## Algoritmo Leading Eigen. Utiliza la idea de modularidad y valores propios del grafo
## para asignar nodos a comunidades. Es eficiente para grafos grandes y complejos.
clust.eigen <- igraph::cluster_leading_eigen(g)$membership



### Evaluando la separación de los clusters - Modularidad ###


## Modularidad es una métrica natural para evaluar la separación entre comunidades/clusters

## La modularidad es una medida que evalúa la calidad de la partición de un grafo en comunidades.
## Cuanto mayor sea la modularidad, mejor será la partición

## La modularidad se define como la diferencia (escalada) entre el peso total observado
## de las aristas entre los nodos en el mismo cluster y el peso total esperado si los
## pesos fueran distribuidos aleatoriamente entre todos los pares de nodos

## Nosotros calcularemos un score de modularidad para cada cluster usando las tasas
## en vez de las diferencias, debido a que las tasas no se ven tan fuertemente
## influenciadas por el tamaño de los clusters !!!


## Obteniendo la métrica de modularidad
library(bluster)

## Matriz de índices de modularidad pairwise
ratio <- pairwiseModularity(g, clust, as.ratio = TRUE) # g es el grafo y clust las (comunidad) en el grafo

dim(ratio)

## Heatmap de los scores de modularidad (ratio)
library(pheatmap)

## Convertimos a logaritmo
## Agregamos +1 para que los 0s no den error por indeterminación
pheatmap(log2(ratio + 1),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100)
)

## Un dataset que contiene clusters bien separados debería contener la mayoría
## del peso total observado en las entradas diagonales, i.e la mayoría de las
## aristas ocurren entre células del mismo cluster


### Ejercicio ###

## Obtener el peso total de la matriz (ratio) en la diagonal
ratio_m_diagonal <- sum(diag(ratio), na.rm = TRUE)

## Obtener el peso total de la matriz (ratio) fuera de la diagonal
ratio_m_total <- sum(ratio, na.rm = TRUE)
ratio_m_total - ratio_m_diagonal




### Otros métodos de clustering ###

## Clustering por k-means ##

# PRO: Rápido
# Se debe especificar el número de clusters de antemano
# Favorece clusters esféricos

## Clustering jerárquico ##

# Produce un dendograma (árbol) representando las células y la similaridad entre
# subpoblaciones a varias resoluciones
# Demasiado lento para correrse en algo más grande que los datasets más pequeños de scRNA-seq



### Evaluando la estabilidad de los clusters ###

## Una propiedad deseable de un cluster dado es que éste sea estable a las perturbaciones
## en los datos de entrada, de esta manera:

# Pequeños cambios al procesamiento no cambiarán el resultado
# Se incrementa la probabilidad de que las conclusiones puedan ser replicadas
# en un estudio independiente

## Valores más altos significan que los clusters se repiten más veces


## Uno puede hacer un proceso de bootstrap para evaluar la estabilidad de un algoritmo
## de clustering en un dataset dado y calcular la coasignación. La coasignación es
## la probabilidad de que células elegidas al azar del cluster X y Y sean asignadas
## al mismo cluster en la réplica del proceso de bootstrap


## Esta función toma un objeto SingleCellExperiment como entrada, realiza una reducción
## de dimensionalidad usando PCA, construye un grafo de vecinos más cercanos basado
## en la similitud de Jaccard, y finalmente asigna clústeres a las células utilizando
## el algoritmo de Louvain.
myClusterFUN <- function(x) {
  g <- buildSNNGraph(x, use.dimred = "PCA", type = "jaccard")
  igraph::cluster_louvain(g)$membership
}

## Aplicamos la función myClusterFUN anteriormente definida a sce.pbmc
originals <- myClusterFUN(sce.pbmc)

## evaluar la estabilidad de los clústeres obtenidos por el algoritmo de Louvain
set.seed(0010010100)
## Hacer el proceso de bootstrap y calcular la probabilidad de coasignación
coassign <- bootstrapStability(sce.pbmc,
                               FUN = myClusterFUN,
                               clusters = originals
)


## Visualizamos co un heatmap
pheatmap(coassign,
         cluster_row = FALSE, cluster_col = FALSE,
         color = rev(viridis::magma(100))
)


## Nuevamente, queremos que los valores más altos estén en la diagonal
## (en este caso los valores de coasignación) !!!

## Probabilidad alta de coasignación indica que X no es estable con respecto a su separación de Y.

## Queremos altas probabilidades de coasignación en la diagonal

## Debes considerar que el bootstraping solo considera el efecto del ruido de muestreo
## e ignora otros factores que pueden afectar la reproducibilidad (como efectos de batch
## o variación entre los donadores)

## Además, una pobre separación puede ser altamente estable



### Subclustering ###

## Mejora la resolucón al repetir el proceso de feature selection y clustering
## dentro de un único cluster

## Se enfoca en los HGVs y PCs que son los más relevantes para un cluster específico

## Se construye un grafo de vecinos más cercanos (SNN) utilizando las coordenadas
## de PCA de las células en sce.pbmc
g.full <- buildSNNGraph(sce.pbmc, use.dimred = "PCA")

## Encontrar los clusters por el algoritmo Walktrap
clust.full <- igraph::cluster_walktrap(g.full)$membership

## Se agrega la asignación de clústeres obtenida a sce.pbmc como clust.full
sce.pbmc$clust.full <- factor(clust.full)

## Visualizar la expresión de los genes "CD3E", "CCR7", "CD69", y "CD44" en diferentes clústeres.
## Las células se agrupan en el eje x (x = "clust.full") y se colorean según su
## asignación de clústeres (colour_by = "clust.full").
plotExpression(sce.pbmc, # plotExpression()
               features = c("CD3E", "CCR7", "CD69", "CD44"),
               x = "clust.full", colour_by = "clust.full"
)

## Conocimiento biológico ##

## CD3E, CCR7, CD69, y CD44 son marcadores de células T de memoria.
## Dentro de las células T de memoria, ¿dónde están las subpoblaciones CD4+ y CD8+?

## El cluster 10 parece tener valores altos para CD3E, CCR7, CD69, y CD44


## Nos quedamos con un subset de columnas (células) correspondientes al cluster 10
memory <- 10

## Nuevo objeto sce.memory que contiene solo las células que pertenecen al clúster número 10
sce.memory <- sce.pbmc[, clust.full == memory]

## modelGeneVar() para calcular las varianzas genéticas modeladas para el clúster 10
dec.memory <- modelGeneVar(sce.memory)

## denoisePCA() para obtener los PCs principales (dentro se obtienen los HVGs)
sce.memory <- denoisePCA(sce.memory,
                         technical = dec.memory, # toma en cuenta la información de variabilidad genética modelada
                         subset.row = getTopHVGs(dec.memory, prop = 0.1) # utiliza solo los HGVs
)


## Se repite el clustering en el subset

## Se construye un grafo de vecinos más cercanos (SNN) utilizando las coordenadas
## de PCA de las células en sce.memory (subset correspondiente al cluster 10)
g.memory <- buildSNNGraph(sce.memory, use.dimred = "PCA")

# Encontrar los clusters por el algoritmo Walktrap
clust.memory <- igraph::cluster_walktrap(g.memory)$membership

## Se agrega la asignación de clústeres obtenida a sce.memory como clust.memory
sce.memory$clust.memory <- factor(clust.memory)

## Visualizar la expresión de los marcadores CD4+ y CD8+ en en los diferentes clústeres  genrados
## para las células del objeto sce.memory  (subset de las células correspondientes al cluster 10)
plotExpression(sce.memory,
features = c("CD8A", "CD4"),
x = "clust.memory"
)

## La expresión de CD4 en este caso es baja, por lo tanto, su cambio es modesto,
## pero la interpretación es clara

## scran::quickSubCluster() ciclará sobre los clusters y realizará el proceso de
## subclustering de acuerdo a una función especificada por el usuario.
## Esto asume que la misma función es apropiada para todos los clusters

## Si tipos celulares o estados celulares se extienden sobre las fronteras de los
## clusters, entonces un subcluster podría representar contaminación de un tipo
## celular en un cluster separado
