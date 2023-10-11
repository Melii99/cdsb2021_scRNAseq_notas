########## Análisis de expresión diferencial ##########

## Paquetes usados ##
library("MouseGastrulationData") ## para descargar datos de ejemplo
library("scater") ## para gráficas y control de calidad
library("scran") ## para selección de genes, clustering, etc
library("batchelor") ## para métodos de correción de batch (lote)
library("patchwork") ## para agrupar gráficas
library("Polychrome") ## para muchos colores
library("bluster") ## métodos de clustering
library("edgeR") ## para expresión diferencial



## scRNA-seq nos puede ayudar a estudiar cambios en composición (cambios en
## proporciones de células) o cambios en niveles de expresión de genes entre
## varias condiciones biológicas


## El primero se llama cambios de abundancia,
# Ejemplo: después de un tratamiento con una droga
# Ejemplo: después de modificaciones genéticas

## Nos permite obtener mayor resolución biológica que experimentos convencionales
## de RNA-seq, sobre todo si podemos asociar cambios en poblaciones celulares a
## manipulaciones experimentales


### Dos categorías de análisis ###


## Análisis de expresión diferencial ##

## Buscamos cambios en niveles de expresión entre condiciones para células del
## mismo tipo que están presentes en todas las condiciones


## Análisis de abundancia diferencial ##

## Buscamos cambios en la composición de los tipos celulares entre condiciones

## Podría ser entre estados celulares en vez de tipos celulares
## Son dos lados de la misma moneda



### Datos de ejemplo ###

## Embriones de ratón quiméricos. Pijuan-Sala, B. et al. A single-cell molecular
## map of mouse gastrulation and early organogenesis. Nature 566, 490–495 (2019).

## Descarga de los datos desde bioconductor y guardar el SCE
library("MouseGastrulationData")
sce.chimera <- WTChimeraData(samples = 5:10)

## Explorar los datos: resumen estadístico para cada columna
sapply(colData(sce.chimera)[, -(1:2)], function(x) {
  x <- if (is.character(x) || is.integer(x)) factor(x) else x
  summary(x)
})

## Encontramos que:

## sample: 6 ratones diferentes
## tomato: inyectados o no con td-Tomato
## pool: lote de secuenciación, cada lote con 1 con y otro sin inyección
## celltype.mappped: 35 tipos de células anotados


## Número de células en nuestras variables principales ##
## (tabla con la distribución de células en función de las combinaciones únicas
## de los valores de sample, pool y tomato)
with(colData(sce.chimera), table(sample, pool, tomato))


## Número de tipos celulares ##
## Cálculo del número de tipos celulares únicos presentes en sce.chimera
length(unique(sce.chimera$celltype.mapped))



### Procesamiento de los datos ###

## Usaremos batchelor porque tenemos muestras de 3 lotes de muestras y
## queremos eliminar diferencias entre los lotes !!!


## Anotacion ##

library("scater")
rownames(sce.chimera) <- uniquifyFeatureNames(
  rowData(sce.chimera)$ENSEMBL, rowData(sce.chimera)$SYMBOL
)


## Control de calidad QC ##

## Filtrar y eliminar las células ya marcadas como "stripped" o "Doublet"
drop <- sce.chimera$celltype.mapped %in% c("stripped", "Doublet")
sce.chimera <- sce.chimera[, !drop]


## Normalizacion ##
sce.chimera <- logNormCounts(sce.chimera)


## Modelado de la varianza ##

## Modelar la varianza yseleccionar los genes con variabilidad biológica distinta de 0 (HVGs)
library("scran")
dec.chimera <- modelGeneVar(sce.chimera, block = sce.chimera$sample)
chosen.hvgs <- dec.chimera$bio > 0


## Merging ##

## corregir batch effects y combinar los datos de diferentes experimentos
library("batchelor")
set.seed(01001001)
merged <- correctExperiments(sce.chimera, #SCE
                             batch = sce.chimera$sample, # batches
                             subset.row = chosen.hvgs, # HVGs
                             PARAM = FastMnnParam( # Orden de los datos
                               merge.order = list( # Orden en que se combinan
                                 list(1, 3, 5), # WT (3 replicates)
                                 list(2, 4, 6) # td-Tomato (3 replicates)
                               )
                             )
)


## Clustering ##

## Creación de un grafo SNN
g <- buildSNNGraph(merged, use.dimred = "corrected")

## Encontrar los clusters aplicando el algoritmo de Louvain
clusters <- igraph::cluster_louvain(g)

## Asignación de etiquetas de colores a las células (merged) con su respectivo clúster
colLabels(merged) <- factor(clusters$membership)


## Reducción de dimensionalidad ##

## Reducción dimensionalidad tSNE
merged <- runTSNE(merged, dimred = "corrected", external_neighbors = TRUE)

## Reducción dimensionalidad UMAP
merged <- runUMAP(merged, dimred = "corrected", external_neighbors = TRUE)



### Exploración de los datos de ejemplo ###


## Exploremos si tenemos clusters con una diferencia grande en el número de celulas
## entre las muestras sin y con inyecciones de td-Tomato

## Exploremos el número de células en cada cluster a lo largo de los
## 3 lotes de secuenciación (batch)


## Clusters vs DE por td-Tomato
table(colLabels(merged), merged$tomato)

## Clusters vs lotes de muestras (batch)
table(colLabels(merged), merged$pool)

## Visualización ##

## Visualizar nuestros clusters que (26 en dimensiones reducidas) de t-SNE
## Queremos que todos los clusters tengan muestras de cada lote de secuenciación (batch).
## Vemos que no parece que haya mucha señal en base a td-Tomato
library("patchwork")
plotTSNE(merged, colour_by = "tomato", text_by = "label") +
  plotTSNE(merged, colour_by = data.frame(pool = factor(merged$pool)))

## Podemos usar facet_wrap() para reducir el over-plotting y ver mejor la información:
##  t-SNE donde las células se agrupan por las categorías de tomato y se visualizan en
## paneles separados para cada categoría
plotTSNE(merged,
         colour_by = "tomato",
         other_fields = c("tomato")
) + facet_wrap(~tomato)

##  t-SNE donde las células se agrupan por las categorías de pool y se visualizan en
## paneles separados para cada categoría
plotTSNE(merged,
         colour_by = data.frame(pool = factor(merged$pool)),
         other_fields = c("pool")
) + facet_wrap(~pool)
