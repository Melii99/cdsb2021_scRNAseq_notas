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



###  Continuación ###


## ¿Algunos de estos genes están asociados con los resultados de clustering?

## P. ej. ¿se encuentra el gen 1 asociado con el clustering?
plotExpression(sce.pbmc,
               features = rownames(sce.pbmc)[1], # features especifica los genes cuya expresión se va a visualizar
               x = "cluster", colour_by = "cluster" # x indica la variable que se utilizará para agrupar las células en el eje x
)


## P. ej. ¿se encuentra el gen 2 asociado con el clustering?
plotExpression(sce.pbmc,
               features = rownames(sce.pbmc)[2],
               x = "cluster", colour_by = "cluster"
)


## P. ej. ¿se encuentra el gen 2512 asociado con el clustering?
plotExpression(sce.pbmc,
               features = rownames(sce.pbmc)[2512],
               x = "cluster", colour_by = "cluster"
)

## P. ej. ¿se encuentra el gen CD3E asociado con el clustering?
plotExpression(sce.pbmc,
               features = "CD3E",
               x = "cluster", colour_by = "cluster"
)


## Ver una gráfica como una forma de encontrar los genes marcadores obviamente
## no nos sirve a gran escala

## Necesitamos un método estadístico para identificar estos genes marcadores

## La prueba t de Welch es una opción obvia para probar las diferencias en
## la expresión entre clústeres !!!



### Prueba t modificada de Welch pareada ###

## Prueba rápida con buenas propiedades estadísticas para un gran número de células

## as comparaciones pareadas proveen un log-fold change para indicar cuáles
## clústerse son distinguidos por cada gen

## ¿Por qué no comparar cada clúster con el promedio de todas las otras células?
## Sensibilidad a la composición poblacional, una subpoblación dominante sola que dirige
## la selección de los marcadores top para cualquier otro clúster



### Ejemplo ilustrativo: CD3E como gen marcador en el dataset PBMC4k 10X ###


## Pruebas pareadas  (combinatoria) ##

## Combinando comparaciones del gen CD3E para el clúster 1
## “Me interesa saber si el gen CD3 está diferencialmente expresado entre el clúster 1 y ..”

# cualquier (any) otro clúster = P = 1.3 x 10-205 (Simes adjusted P-value)
# todos (all) los otros clústeres = P = 0.11 (Berger’s intersection-union test)
# algunos (some) de los otros clústeres = P = 2.0 x 10-44 (mediana u otro cuantil, Holm-adjusted P-values)


## Extendiendo a todos los genes - Aplicación estándar ##

## Para cada clúster, usar pruebas t de Welch para identificar los genes que están
## diferencialmente expresados entre éste y cualquier (any) otro clúster

## Usadas para todos los genes
#scran::pairwiseTTests()
#scran::combineMarkers()

## Usando findMarkers() encontrar marcadores diferenciales entre grupos en datos de células individuales
library(scran)

markers.pbmc <- findMarkers(sce.pbmc,
                            groups = sce.pbmc$cluster, # EWncontrar marcadores agrupando por clusters
                            test.type = "t", pval.type = "any" # Prueba t, comparación con todos los clusters
)

## La lista resultante contiene dataframes correspondientes a cada cluster
## Cada dataframe (genes marcadores) ya se encuentran ordenados por p-value (más significativos)
## En caso de querer ordenar por logfold change (genes más sobreexpresados o subexpresados)
## es posible pero puede que no sean tan significativos


## Explorando los resultados ##

## Seleccion del clúster 9 para explorar los genes marcadores asociados con este clúster
chosen <- "9"

## Se almacenan los genes marcadores del clúster seleccionado en la variable interesting
interesting <- markers.pbmc[[chosen]]

## plotExpression() para visualizar la expresión de los cuatro primeros genes marcadores
plotExpression(sce.pbmc, rownames(interesting)[1:4],
               x = "cluster", colour_by = "cluster"
)

## Explorar los genes marcadores obtenidos con un heatmap
best.set <- interesting[interesting$Top <= 6, ]
logFCs <- as.matrix(best.set[, -(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))
library(pheatmap)
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))

## Usamos el campo Top para identificar un conjunto de genes que distinguen el
## clúster 9 de cualquier otro clúster !!!



### Log-fold change ###


## Sin espeficiar el log-fold change ##

## Para cada clúster, usa pruebas t de Welch para identificar los genes que están
## sobreexpresados entre éste y cualquier otro clúster

## Identificar los genes que tienen una expresión significativamente mayor que otros clústeres
markers.pbmc.up <- findMarkers(sce.pbmc,
                               groups = sce.pbmc$cluster,
                               test.type = "t", direction = "up", pval.type = "any"
)
## Almacenar los genes upregulated
interesting.up <- markers.pbmc.up[[chosen]]


## Especificando el lfc ##

## Para cada clúster, usa pruebas t de Welch para identificar los genes que están
## sobreexpresados con un log-fold change (lfc) o al menos 1 entre éste y cualquier otro clúster

## Identificar los genes que tienen una expresión significativamente mayor que otros clústeres
## Pero ahora, con lfc = 1 se establece un punto de corte (reduciendo la lista)
markers.pbmc.up2 <- findMarkers(sce.pbmc,
                                groups = sce.pbmc$cluster,
                                test.type = "t", direction = "up", lfc = 1, pval.type = "any"
)

interesting.up2 <- markers.pbmc.up2[[chosen]]

## La prueba t también nos permite especificar un log-fold change diferente de
## cero como la hipótesis nula

## Es más riguroso que simplemente filtrar por log-fold change TREAT


## Heatmap ##

best.set <- interesting.up2[interesting.up2$Top <= 5, ]
logFCs <- as.matrix(best.set[, -(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))


## Los promedios están más centrados en un conjunto de genes marcadores candidatos
## que están sobreexpresados en el clúster 9

## El incremento del rigor no se da sin costo

## Si el lfc es muy grande podría descartar genes útiles (punto de corte muy estricto)

## E.g., un gen sobreexpresado en una proporción pequeña de células en un clúster sigue
## siendo un marcador efectivo si el foco está en la especificidad más que en la sensibilidad



### Encontrando marcadores específicos de clústeres ###


## Por defecto, scran::findMarkers() dará un alto rango a genes que están DE en
## cualquier comparación pareada

## Quiero genes que son específicos de cada clúster

## Tú quieres genes que son DE en todas las comparaciones pareadas
## Para cada clúster, usa pruebas t de Welch para identificar genes que están
## sobreexpresados entre éste y todos los otros clústeres


## Identificar genes que tienen una expresión significativamente mayor en
## comparación con otros clústeres, y usando pval.type = "all" se considerarán
## todos los valores p, independientemente de si están ajustados o no
## Es decir, se identifican genes solo upregulated para ese cluster
markers.pbmc.up3 <- findMarkers(sce.pbmc,
                                groups = sce.pbmc$cluster,
                                direction = "up", pval.type = "all" #
)

##
interesting.up3 <- markers.pbmc.up3[[chosen]]

## CONTRAS:

## Los genes pueden aparecer con sobreexpresión para diferentes tipos celulares
## (Ej. poblaciones DN(CD4-/CD8-) no hay expresión en CD4 NI EN CD8 ) por lo que
## con sólo uno no siempre se pueden distinguir (a menos que las poblaciones sean
## muy diferentes y esten bien definidas)



## findMarkers con pval.type some ##


## Útil cuando pval.type="all" es muy estricto todavía pval.type="any" es muy generoso

## Aplica la corrección Holm-Bonferroni a los P-values y toma el mejor valor de
## en medio como el P-value combinado

## Perderás algunas garantías ofrecidas por los otros métodos !!!


## Para cada clúster, usa pruebas t de Welch para identificar los genes que están
## sobreexpresados entre éste y algunos de los otros clústers
markers.pbmc.up4 <- findMarkers(sce.pbmc,
                                groups = sce.pbmc$cluster,
                                direction = "up", pval.type = "some"
)

interesting.up4 <- markers.pbmc.up4[[chosen]]
