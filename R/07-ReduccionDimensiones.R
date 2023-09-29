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



### Juntando todo ###

set.seed(100)
## Generar y almacenar el set completo de componentes principales  con los HVGs
sce.zeisel <- runPCA(sce.zeisel, subset_row = top.zeisel)

## También pueden seleccionarse un subset de componentes principales usando
## el metodo del codo

# library(PCAtools)
# percent.var <- attr(reducedDim(sce.zeisel), "percentVar")
# chosen.elbow <- PCAtools::findElbowPoint(percent.var)
reducedDim(sce.zeisel, "PCA_elbow") <- reducedDim(
  sce.zeisel, "PCA"
)[, 1:chosen.elbow]

## O seleccionar un subset de componentes principales usando métodos basados
## en la estructura poblacional

# choices <- getClusteredPCs(reducedDim(sce.zeisel))
# chosen.clusters <- metadata(choices)$chosen
reducedDim(sce.zeisel, "PCA_clusters") <- reducedDim(
  sce.zeisel, "PCA"
)[, 1:chosen.clusters]


### Ejercicio: PCA de sce.pbmc ###

## Realiza un PCA para los datos sce.pbmc.
## Elige el número de PCs a conservar utilizando el método del codo
## Elige el número de PCs a conservar utilizando la estructura de la población
## Agrega esta información al objeto sce.pbmc


## PCA de los topGenes obtenidos para sce.pbmc
library(scran)
library(scater)

set.seed(100)
sce.pbmc <- runPCA(sce.pbmc,
                     subset_row = top.pbmc
)

## Encontrar el numero óptimo de PCs con el método del codo
library(PCAtools)

percent.var2 <- attr(reducedDim(sce.pbmc), "percentVar") # reducedDim()
chosen.elbow2 <- PCAtools::findElbowPoint(percent.var2) # findElbowPoint()

## Encontrar el numero óptimo de PCs utilizando la estructura de la población
choices2 <- getClusteredPCs(reducedDim(sce.pbmc)) # getClusteredPCs()
chosen.clusters2 <- metadata(choices2)$chosen

## Agregar información del numero óptimo de PCs con el método del codo a sce.pbmc
library(PCAtools)

reducedDim(sce.pbmc, "PCA_elbow") <- reducedDim(
  sce.pbmc, "PCA"
)[, 1:chosen.elbow2]

## Ver la información agregada
head(reducedDim(sce.pbmc, "PCA_elbow"))

## Agregar información del numero óptimo de PCs utilizando la estructura de la población a sce.pbmc
reducedDim(sce.pbmc, "PCA_clusters") <- reducedDim(
  sce.pbmc, "PCA"
)[, 1:chosen.clusters2]

## Ver la información agregada
head(reducedDim(sce.pbmc, "PCA_clusters"))



### Usando el ruido técnico ###

## Otra alternativa es conservar todos los PCs hasta que el % de variación explicado
## alcance algun límite (por ejemplo, basado en la estimación de la variación técnica)

## denoisePCA() automáticamente selecciona el número de PCs

## Por default, denoisePCA() realiza algunas simulaciones, por lo tanto necesitamos
## ‘configurar la semilla’ para obtener resultados reproducibles

library(scran)

set.seed(111001001)
denoised.pbmc <- denoisePCA(sce.pbmc,
                            technical = dec.pbmc, subset.row = top.pbmc
)

## Dimensiones
dim(reducedDim(denoised.pbmc, "PCA"))

## La dimensionalidad del output es el límite inferior para el número de PCs requeridos
## para explicar toda la variación biológica. Lo que significa que cualquier número
## menor de PCs definitivamente descartará algún aspecto de la señal biológica

## Esto no grantiza que los PCs retenidos capturen toda la señal biológica

## Esta técnica usualmente retiene más PCs que el método del punto del codo

## scran::denoisePCA() internamente limita el numero de PCs, por default 5-50,
## para evitar la selección de excesivamente pocos PCs (cuando el ruido técnico
##es alto relativo al ruido biológico) o exceso de PCs (cuando el ruido técnico es demasiado bajo




## ¿Qué pasa si calculamos la relación media-varianza con la función modelGeneVar
## para el dataset sce.pbmc? ##

dec.pbmc2 <- modelGeneVar(sce.pbmc)
denoised.pbmc2 <- denoisePCA(sce.pbmc,
                             technical = dec.pbmc2, subset.row = top.pbmc
)
dim(reducedDim(denoised.pbmc2))

## scran::denoisePCA() tiende a funcionar mejor cuando la relación media-varianza
## refleja el ruiudo técnico verdadero, i.e estimado por scran::modelGeneVarByPoisson()
## o scran::modelGeneVarWithSpikes() en vez de scran::modelGeneVar()

## El dataset PBMC está cerca de este límite inferior: el ruido técnico es alto
## relativo al ruido biológico



## ¿Qué pasa si calculamos el número de PCs usando el ruido técnico para el dataset sce.pbmc? ##

set.seed(001001001)
denoised.zeisel <- denoisePCA(sce.zeisel,
                              technical = dec.zeisel, subset.row = top.zeisel
)
dim(reducedDim(denoised.zeisel))

## Los datos de cerebro de Zeisel están cerca de este límite superior: el ruido técnico es demasiado bajo



### Reducción de dimensionalidad para visualización ###

## Clustering y otros algoritmos operaran fácilmente sobre 10-50 (a lo más) PCs,
## pero ese número es aún demasiado para la visualización

## Por lo tanto, necesitamos estrategias adicionales para la reducción de
## dimensionalidad si queremos visualizar los datos



### Visualizando con PCA ###

## PCA es una técnica lineal, por lo tanto, no es eficiente para comprimir
## diferencias en más de 2 dimensiones en los primeros 2 PCs

## Visualizar PC1 y PC2 del PCA de zeisel
plotReducedDim(sce.zeisel, dimred = "PCA")

## Visualizar PC1 y PC2 del PCA de zeisel coloreando los clusters por level1class
## (tipo celular)
plotReducedDim(sce.zeisel,
               dimred = "PCA",
               colour_by = "level1class"
)


## Retos y resumen de la visualización con PCA ##

## PROS:
## PCA es predecible y no introducirá estructura aritficial en los datos.
## Es determínistico y robusto a cambios pequeños en los valores de entrada.

## CONTRAS:
## Usualmente no la visualización no es suficiente para visualizar la naturaleza
## compleja de los datos de scRNA-seq

plotReducedDim(sce.zeisel,
               dimred = "PCA",
               ncomponents = 4, colour_by = "level1class"
)


### Visualización con t-SNE ###

## t-stochastic neighbour embedding (t-SNE) es la visualización por excelencia
## de datos de scRNA-seq. Intenta encontrar una representación (no-lineal) de los
## datos usando pocas dimensiones que preserve las distancias entre cada punto y
## sus vecinos en el espacio multi-dimensional

set.seed(00101001101)
sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA")
plotReducedDim(sce.zeisel, dimred = "TSNE", colour_by = "level1class")

##  Retos de la visualización con t-SNE ##

set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel,
                      dimred = "PCA",
                      perplexity = 30 # número de vecinos cercanos aconsiderar
)

plotReducedDim(sce.zeisel,
               dimred = "TSNE",
               colour_by = "level1class"
)



## Preguntas ##

## ¿Qué pasa si vuelves a correr runTSNE() sin especificar la semilla?
## R. el resultado es aleatorio
## ¿Qué pasa si especificas la semilla pero cambias el valor del parámetro perplexity?
## R. el numero de vecinos a considerar cambia y por ende también el resultado

## Ajustando la perplejidad de tSNE

set.seed(100)

sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA", perplexity = 5)
p1 <- plotReducedDim(sce.zeisel, dimred = "TSNE", colour_by = "level1class")

sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA", perplexity = 20)
p2 <- plotReducedDim(sce.zeisel, dimred = "TSNE", colour_by = "level1class")

sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA", perplexity = 80)
p3 <- plotReducedDim(sce.zeisel, dimred = "TSNE", colour_by = "level1class")

library("patchwork")
p1 + p2 + p3


## Una perplejidad alta permite mejor representación de estructuras globales,
## especialmente para datos que tienen múltiples escalas y estructuras a diferentes distancias.

## Selección de parámetros para t-SNE con cierta profundidad:

## Algunos componentes aleatorios y la selección de parámetros cambiarán la visualización
## La interpretación puede ser engañada por el tamaño y posición de los clusters
## t-SNE infla clusters densos y comprime clusters escasos
## t-SNE no está obligado a preservar las localizaciones relativas de clusters no-vecinos (no puedes interpretar distancias no locales)

## Aún así: t-SNE es una herramienta probada para visualización general de datos de
## scRNA-seq y sigue siendo muy popular



### Visualización con UMAP ###

## Uniform manifold approximation and project (UMAP) es una alternativa a t-SNE

## Así como t-SNE, UMAP intenta encontrar una representación (no lineal) de pocas d
## imensiones de los datos que preserve las distancias entre cada puntos y sus vecinos
## en el espacio multi-dimensional

## t-SNE y UMAP están basados en diferentes teorías matemáticas

set.seed(100)
sce.zeisel <- runUMAP(sce.zeisel,
                      dimred = "PCA",
                      n_neighbors = 15
)
plotReducedDim(sce.zeisel,
               dimred = "UMAP",
               colour_by = "level1class"
)


## Comparado con t-SNE: ##

## UMAP tiende a encontrar clusters visualmente más compactos
## Intenta preservar más de la estructura global que t-SNE
## Tiende a ser más rápido que t-SNE, lo cual puede ser importante para datasets grandes.
## La diferencia desaparece cuando se aplican sobre los primeros PCs


## Retos de la visualización con UMAP ##

## gual que para t-SNE, es necesario configurar una semilla y diferentes valores
## para los parámetros cambiarón la visualización

## Si el valor para los parámetros n_neighbors o min_dist es demasiado bajo entonces
## el ruido aleatorio se interpretará como estructura de alta-resolución,
## si son demasiado altos entonces se perderá la estructura fina

## TIP: Trata un rango de valores para cada parámetro para asegurarte de que no
## comprometen ninguna de las conclusiones derivadas de la gráfica UMAP o t-SNE


## Preguntas ##

## ¿Qué pasa si vuelves a correr runUMAP() sin especificar la semilla?
## R. Al igual que tSNE los resultados son aleartorios y cambila visualización
## ¿Qué pasa si especificas la semilla pero cambias el valor del parámetro n_neighbors?
## R. R. el numero de vecinos a considerar cambia y por ende también el resultado



### Interpretación de la visualización ###

## Reducción de dimensionalidad para la visualización de los datos necesariamente #
## involucra descartar información y distorsionar las distancias entre las células

## No sobre interpretes las gráficas bonitas !!!


### Resumen y recomendaciones ###

## Las gráficas de t-SNE y UMAP son herramientas diagnóstico importantes, por ejemplo:
## para checar si dos clusters son realmente subclusters vecinos o si un cluster puede
## ser dividido en más de un cluster

## Es debatible cual visualización, t-SNE o UMAP, es más útil o estáticamente agradable.

## Está bien elegir aquella que funcione mejor para tu análisis (tomando en cuenta
## que tratarás la gráfica únicamente como una herramienta de visualización/diagnóstico
## y que no llegarás a ninguna conclusión fuerte basado únicamente en la gráfica )
