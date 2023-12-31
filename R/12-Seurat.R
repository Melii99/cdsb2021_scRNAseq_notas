########## Introducción a Seurat ##########


### Una perspectiva diferente ###

## Seurat es un paquete R diseñado para control de calidad, análisis y exploración
## de datos de secuencia de ARN de una sola célula. Seurat tiene como objetivo
## permitir a los usuarios identificar e interpretar fuentes de heterogeneidad
## a partir de mediciones transcriptómicas unicelulares e integrar diversos
## tipos de datos unicelulares.

## Seurat es desarrollado y mantenido por el laboratorio de Satija y se publica
## bajo la Licencia Pública GNU (GPL 3.0).


## En este tutorial se ve como procesar los datos de scRNAseq con un nuevo paquete.
## Los pasos a realizar son en esencia los mismos que ya revisamos con el tutorial
## de la OSCA de RStudio.

## El paquete mas adecuado y que deberás utilizar dependerá mayoritariamente de
## tus datos y el procesamiento que se adecúe a estos.


## Cargar paquetes de R
library("BiocFileCache") ## para descargar datos
library("dplyr") ## para filtar datos
library("Seurat") ## paquete principal de este capítulo
library("patchwork") ## para graficar imágenes juntas



### Dataset de ejemplo ###

## Peripheral Blood Mononuclear Cells (PBMC) disponibles gratuitamente de 10X Genomics.
## Son en total 2,700 céluas únicas secuenciadas con Illumina NextSeq 500.


# Usemos datos de pbmc3k tal y como lo hacen en
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# pero con nuestro propio código

## DEscarga de los datos
bfc <- BiocFileCache()
raw.path <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples",
  "cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
))

untar(raw.path, exdir = file.path(tempdir(), "pbmc3k"))

fname <- file.path(tempdir(), "pbmc3k/filtered_gene_bc_matrices/hg19")

## Leer el dataset PBMC  en R
pbmc.data <- Read10X(data.dir = fname)


## Inicializar el objeto Seurat con los datos crudos ()raw data, no normalizada ##
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

## Imprimir el objeto creado
pbmc


## Estructura del objeto
str(pbmc)

## Ahora, ¿Cómo accedemos a cada slot?

## La forma / dimensiones del objeto se pueden encontrar usando las funciones dim(),
## ncol() y nrow(); Los nombres de celda y característica se pueden encontrar usando
## las funciones colnames() y rownames(), respectivamente, o la función dimnames().
## Se puede obtener un vector de nombres de objetos Assay, DimReduc y Graph contenidos
## en un objeto Seurat mediante el uso de nombres.

dim(pbmc) # Dimensiones

head(rownames(pbmc)) # Primeros rownames (genes)

head(colnames(pbmc)) # Primeros colnames ()


## Se puede obtener un vector de nombres de objetos Assay, DimReduc y Graph
## contenidos en un objeto Seurat mediante el uso de nombres.

## La extracción de objetos específicos de Assay, DimReduc o Graph se puede realizar
## con el operador doble [[ ]] extract. La adición de nuevos objetos a un objeto de
## Seurat también se hace con el operador doble [[ ]] extract; Seurat averiguará a
## qué parte del objeto Seurat pertenece un nuevo objeto asociado.

names(pbmc) # Nombres disponibles en el objeto

pbmc[["RNA"]] # Acceder a RNA

## El acceso a los datos de un objeto Seurat se realiza con la función GetAssayData().
## La adición de datos de expresión a counts, data, o scale.data se puede hacer con
## SetAssayData(). Los datos nuevos deben tener las mismas celdas en el mismo orden
## que los datos de la expresión actual. Los datos agregados a los recuentos o datos
## deben tener las mismas características que los datos de la expresión actual.

GetAssayData(object = pbmc, slot = "data")[1:3, 1:3]

## Metadata de las Células.

## Se puede acceder a los metadatos a nivel de celda con el operador de extracción
## [[ ]] extract o usando $sigil. Extraer con $sigil significa que solo se puede
## extraer un bit de metadatos a la vez, aunque se ha habilitado el autocompletado
## de pestañas, lo que lo hace ideal para uso interactivo. La adición de metadatos
## a nivel de celda se puede configurar usando el operador de extracción único
## [[ ]] también, o usando AddMetaData.

head(pbmc@meta.data)

head(pbmc[[c("nCount_RNA", "nFeature_RNA")]])

## drop = TRUE convertirá los metadatos en un vector de nombres, donde cada entrada
## estará nombrada según la célula a la que corresponde.
head(pbmc[["nCount_RNA", drop = TRUE]])


## La clase Assay almacena datos de una sola celda.

## Para los experimentos típicos de scRNA-seq, un objeto Seurat tendrá un único
## ensayo (“RNA”). Este ensayo también almacenará múltiples ‘transformaciones’
## de los datos, incluidos recuentos sin procesar (ranura @counts), datos normalizados
## (ranura @data) y datos escalados para la reducción dimensional (ranura @scale.data).

## Para experimentos más complejos, un objeto podría contener múltiples ensayos.
## Estos podrían incluir tipos de datos multimodales (etiquetas derivadas de
## anticuerpos CITE-seq, ADT) o mediciones imputadas / corregidas por lotes.
## Cada uno de esos ensayos tiene la opción de almacenar también las mismas
## transformaciones de datos.

## ¿Cómo se ven los datos en una matriz de recuento?

## Examinemos algunos genes en las primeras treinta células. Los valores en la matriz
## representan ceros (no se detectan moléculas). Dado que la mayoría de los valores
## en una matriz scRNA-seq son 0, Seurat utiliza una representación de matriz dispersa
## (sparse matrix) siempre que sea posible. Esto da como resultado un ahorro
## significativo de memoria y velocidad.

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# tamaño en memoria de la matriz pbmc.data después de ser convertida a una matriz densa (dense matrix).
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

# tamaño en memoria de la matriz pbmc.data antes de ser convertida a una matriz densa
sparse.size <- object.size(pbmc.data)
sparse.size

dense.size / sparse.size # ¡EN ESTE CASO UNA MATRIZ NO DISPERSA OCUPA ~ 24 VECES MAS ESPACIO!


### Control de calidad QC ###

## Algunas métricas de control de calidad comúnmente utilizadas por la comunidad incluyen:

# El número de genes únicos detectados en cada célula.

# Las células de baja calidad o las gotitas vacías suelen tener muy pocos genes.

# Los dobletes o multipletes celulares pueden exhibir un recuento de genes aberrantemente alto

# De manera similar, el número total de moléculas detectadas dentro de una célula
# (se correlaciona fuertemente con genes únicos).

# El porcentaje de lecturas que se asignan al genoma mitocondrial.

# Las células de baja calidad / moribundas a menudo exhiben una extensa contaminación mitocondrial


## Calculamos métricas de control de calidad mitocondrial con la función PercentageFeatureSet(),
## que calcula el porcentaje de recuentos que se originan a partir de un conjunto de características.

## El operador [[ puede agregar columnas a los metadatos del objeto. Este es un gran
## lugar para almacenar estadísticas de control de calidad. Entonces calculamos y
## añadimos la cantidad de lecturas que corresponden al genoma mitocondrial.

## Porcentaje de moléculas de mtADN en cada célula con PercentageFeatureSet() y
## almacenarlas en 'percent.mt'
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## Visualizamos las métricas de control de calidad mencionadas anteriormente como
## un diagrama de violín. Ademas vemos la correlación entre el numero de moléculas
## de RNA detectadas en cada célula con el número de genes únicos y con el
## porcentaje de lecturas que corresponden a mtADN.

## Diagrama de violín del número de genes únicos detectados (nFeature_RNA),
## el número total de moléculas de RNA (nCount_RNA), y el porcentaje de moléculas
## de RNA mitocondrial (percent.mt).
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## scatter plot para visualizar la relación entre el número total de moléculas de RNA
## (nCount_RNA) y el porcentaje de mtADN (percent.mt).
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

## scatter plot para para visualizar la relación entre el número total de moléculas
## de RNA (nCount_RNA) y el número de genes únicos detectados (nFeature_RNA).
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## Combinar los dos gráficos de dispersión
plot1 + plot2


## Finalmente filtramos aquellas células que se salen de los estándares de cada uno de los parámetros.

## Conservar las células con más de 200 genes únicos detectados, menos de 2500 genes
## únicos detectados y un porcentaje MT inferior al 5%.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


## ¿Dónde se almacenan la métricas de QC en Seurat?
## Están almacenadas en la seccion de @meta.data del objeto Seurat.
head(pbmc@meta.data, 5)



### Normalización ###

## De forma predeterminada, se emplea un método de normalización de escala global
## "LogNormalize” que normaliza las medidas de expresión de características para
## cada célula por la expresión total, multiplica esto por un factor de escala
## (10.000 por defecto) y transforma el resultado en logaritmos. Los valores
## normalizados se almacenan en pbmc [["RNA"]]@data.

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)



### Detección de genes (caractersticas) altamente variables ###


## A continuación, calculamos un subconjunto de características que exhiben una
## alta variación de célula a célula en el conjunto de datos (es decir, están altamente
## expresadas en algunas células y poco expresadas en otras). El equipo de Seurat
## y otros equipos han descubierto que centrarse en estos genes en el análisis
## posterior ayuda a resaltar la señal biológica en conjuntos de datos unicelulares.

## El procedimiento en Seurat mejora a comparación de las versiones anteriores
## al modelar directamente la relación de varianza media inherente a los datos de
## una sola célula, y se implementa en la función FindVariableFeatures().
## De forma predeterminada, se devuelven 2000 características por conjunto de datos
## (aunque se puede modificar). Estos se utilizarán en análisis posteriores, como PCA.

## Encomtramos HVGs con FindVariableFeatures()
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## Identificar los 10 genes más altamente variables
top10 <- head(VariableFeatures(pbmc), 10)

## Top 10 HVGs
top10

## Graficar los HVGs con y sin labels
library(Seurat)
VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = top10, repel = TRUE)



### Escalar los datos ###


## A continuación, aplicamos una transformación lineal (“escalado”) que es un paso
## de preprocesamiento estándar antes de las técnicas de reducción dimensional como PCA.


## La función ScaleData():

# Cambia la expresión de cada gen, de modo que la expresión media en las células sea 0

# Escala la expresión de cada gen, de modo que la varianza entre las células sea 1

# Este paso otorga el mismo peso en los análisis posteriores, de modo que los genes
# altamente expresados no dominen

# Los resultados de esto se almacenan en pbmc [["RNA"]]@scale.data.

## Escalar todos los genes con ScaleData()
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)



### Reducción dimensional lineal - PCA ###


## A continuación, realizamos PCA sobre los datos escalados. De forma predeterminada,
## solo las características variables determinadas previamente se utilizan como entrada,
## pero se pueden definir mediante el argumento de características si desea elegir
## un subconjunto diferente.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

## Seurat proporciona varias formas útiles de visualizar tanto las células como
## las características que definen el PCA, incluidas VizDimReduction(), DimPlot()
## y DimHeatmap().


## Es posible examinar y visualizar los resultados de PCA de diferentes formas:

## resultados del análisis de PCA para las dimensiones 1 a 5 y muestra las características
## (genes) más importantes asociadas con cada dimensión.
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

## visualización de los pesos (cargas) de las características para las dos
## primeras dimensiones del PCA.
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

## gráfico de dispersión para visualizar las células en las dos primeras dimensiones del PCA
DimPlot(pbmc, reduction = "pca")


## En particular, DimHeatmap() permite una fácil exploración de las fuentes
## primarias de heterogeneidad en un conjunto de datos y puede ser útil cuando se
## intenta decidir qué PC incluir para análisis posteriores posteriores. Tanto las
## células como las características se ordenan de acuerdo con sus puntajes de PCA.
## Establecer cells en un número traza las células “extremas” en ambos extremos del
## espectro, lo que acelera drásticamente el trazado de grandes conjuntos de datos.
## Aunque claramente es un análisis supervisado, consideramos que esta es una
## herramienta valiosa para explorar conjuntos de características correlacionadas.


## Heatmap  para la primera dimensión  del PCA. Muestra la expresión génica de las 500
## células con las proyecciones más altas (en la primera dimensión del PCA).
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

## Heatmaps (15)  para las primera s 15 dimensión  del PCA. Muestra la expresión génica
## de las 500  células con las proyecciones más altas en cada dimensión
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)



### Determinar la dimensionalidad del conjunto de datos ###


## Para superar el extenso ruido técnico en cualquier característica única para
## los datos de scRNA-seq, Seurat agrupa las células en función de sus puntuaciones
## de PCA, y cada PC representa esencialmente una “metafunción” que combina información
## en un conjunto de características correlacionadas. Por lo tanto, los componentes
## principales principales representan una compresión sólida del conjunto de datos.
## Sin embargo, ¿cuántos componentes deberíamos elegir incluir? 10? 20? 100?


# En Macosko et al, implementamos una prueba de remuestreo inspirada en el procedimiento
# JackStraw. Permutamos aleatoriamente un subconjunto de los datos (1% por defecto) y
# volvemos a ejecutar PCA, construyendo una “distribución nula” de puntuaciones de
# características, y repetimos este procedimiento. Identificamos PC “importantes”
# como aquellas que tienen un gran enriquecimiento de características de bajo valor p.


# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

## La función JackStrawPlot() proporciona una herramienta de visualización para
## comparar la distribución de los p-values para cada PC con una distribución uniforme
## (línea discontinua). Las PC “significativas” mostrarán un gran enriquecimiento de
## funciones con valores p bajos (curva sólida por encima de la línea discontinua).
## En este caso, parece que hay una fuerte caída en la importancia después de los
## primeros 10-12 PCs.

## Gráfico que compara la distribución de los valores p para cada PC (1:15) con una
## distribución uniforme (línea discontinua)
JackStrawPlot(pbmc, dims = 1:15)

## Un método heurístico alternativo genera un “diagrama de codo (Elbow Plot)”:
## una clasificación de componentes principales basada en el porcentaje de varianza
## xplicada por cada uno (función ElbowPlot()). En este ejemplo, podemos observar un
## “codo” alrededor de PC9-10, lo que sugiere que la mayor parte de la señal verdadera
## se captura en las primeras 10 PC.

## Visualización del método elbow para notar la "caida" en el porcentaje de varianza
## explicada por cada PC
ElbowPlot(pbmc)



### Clustering ###


## Construcción de un grafo SNN usando los primeros 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

## Encontrar clusters usando el algoritmo de Louvain
pbmc <- FindClusters(pbmc, resolution = 0.5)

## resolution es un parámetro que controla la resolución del agrupamiento:
## Un valor más alto de resolution produce un número mayor de grupos, mientras
## que un valor más bajo produce un número menor de grupos más grandes.



### Reducción dimensional no lineal (UMAP/tSNE) ###


## Seurat ofrece varias técnicas de reducción dimensional no lineal, como tSNE y
## UMAP, para visualizar y explorar estos conjuntos de datos. El objetivo de estos
## algoritmos es aprender la variedad subyacente de los datos para colocar células
## similares juntas en un espacio de baja dimensión. Las células dentro de los grupos
## basados en gráficos determinados anteriormente deben ubicarse conjuntamente en
## estos gráficos de reducción de dimensión. Como entrada para UMAP y tSNE, se sugiere
## usar las mismas PC como entrada para el análisis de agrupamiento.


## Generar UMAP con los promeros 10 PCs
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

## Visualizar el UMAP
# Note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


## Se puede guardar el objeto en este punto para que se pueda volver a cargar
## fácilmente sin tener que volver a ejecutar los pasos computacionalmente intensivos
## realizados anteriormente o compartir fácilmente con los colaboradores.

#if (interactive()) {
#  saveRDS(pbmc, file = "pbmc_tutorial.rds")
#}



### Caracteristicas diferencialmente expresadas (biomarcadores de los clusters) ###


## Seurat puede ayudar a encontrar marcadores que definan clústeres mediante
## expresión diferencial. De forma predeterminada, identifica marcadores positivos
## y negativos de un solo grupo (especificado en ident.1), en comparación con todas
## las demás células. FindAllMarkers() automatiza este proceso para todos los clústeres,
## pero también se pueden comparar grupos de clústeres entre sí o contra todas las células.

## El argumento min.pct requiere que se detecte una característica en un porcentaje
## mínimo en cualquiera de los dos grupos de células, y el argumento thresh.test
## requiere que una característica se exprese diferencialmente (en promedio) en
## alguna cantidad entre los dos grupos. Puede establecer ambos en 0, pero con un
## aumento dramático en el tiempo, ya que esto probará una gran cantidad de
## características que probablemente no sean altamente discriminatorias.

## ¿Demasiado lento?

## Como otra opción para acelerar estos cálculos, se puede configurar el número
## máximo de células por identificador. Esto reducirá la resolución de cada clase
## de identidad para que no tenga más células que las que se establezcan.
## Si bien generalmente habrá una pérdida de potencia, los aumentos de velocidad
## pueden ser significativos y es probable que las características expresadas de
## manera más diferencial aún se eleven a la cima.


## Encontrar todos los marcadores para el cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

## Encontrar todos los marcadoresque distinguen al cluster 5 de los clusters 0 y 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

## Encontrar genes marcadores (sobreexpresión) para cada cluster comparados con todas las células restantes

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# only.pos = TRUE - se reportan solo los positivos (los genes que están sobreexpresados)
# min.pct = 0.25 - porcentaje mínimo de células
# logfc.threshold = 0.25 - umbral de cambio en el logfold change
# (su expresión es al menos un 25% mayor en un grupo en comparación con otro grupo.)


## utilizar el paquete dplyr en R para manipular y analizar los resultados de marcadores
## (pbmc.markers) obtenidos previamente

## Selección de los dos mejores marcadores para cada cluster en función de su avg_log2FC
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)


## Seurat tiene varias pruebas de expresión diferencial que se pueden configurar
## con el parámetro test.use. Por ejemplo, la prueba ROC devuelve el “poder de clasificación”
## para cualquier marcador individual (que varía de 0 - aleatorio a 1 - perfecto)

## Encontrar marcadores específicos para el cluster 0
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# ident.1 = 0 - 0 es el identificador del cluster que estás seleccionando.
# logfc.threshold = 0.25 - umbral de fold change mínimo para considerar un gen como marcador
# test.use = "roc" - Realizar prueba de curva ROC para evaluar la capacidad de discriminación de los genes
# only.pos = TRUE - Solo se consideran los genes que tienen una sobreexpresión como marcadores


## Se incluyen varias herramientas para visualizar la expresión de los marcadores.
## VlnPlot() (muestra distribuciones de probabilidad de expresión entre clústeres)
## y FeaturePlot() (visualiza la expresión de características en un gráfico tSNE o
## PCA) son nuestras visualizaciones más utilizadas. También sugerimos explorar RidgePlot(),
## CellScatter() y DotPlot() como métodos adicionales para ver su conjunto de datos.


## Gráfico de violín para las expresiones de los genes MS4A1 y CD79A
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))


## Gráfico de violín para los genes NKG7 y PF4 utilizando las cuentas crudas
## (raw counts) en lugar de los valores de expresión transformados.
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

## Gráfico de puntos para la expresión de los genes indicados (MS4A1, GNLY, CD3E, etc.)
## en diferentes células (permite visualizar la distribución de expresión en diferentes clusters)
FeaturePlot(pbmc, features = c(
  "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
  "CD8A"
))


# DoHeatmap() generates an expression heatmap for given cells and features. In this case,
# we are plotting the top 10 markers (or all markers if less than 10) for each cluster.

##Heatmap de expresión para los 10 genes principales (marcadores) en cada cluster
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()



### Asignar la identidad del tipo celular a los clusters ###


## Podemos usar marcadores canónicos para hacer coincidir fácilmente la agrupación
## imparcial con los tipos de células conocidos.

## asignando nuevos nombres a los clusters en tus datos de scRNA-seq y visualizarlos
## en un gráfico de dispersión

## Vector con los nuevos nombres de los clusters en el mismo orden que los niveles del objeto pbmc.
new.cluster.ids <- c(
  "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
  "NK", "DC", "Platelet"
)

names(new.cluster.ids) <- levels(pbmc) # mismo orden

## usando la función RenameIdents(), asignas estos nuevos nombres a los clusters en el objeto pbmc.
pbmc <- RenameIdents(pbmc, new.cluster.ids)

## Scaterplot con nuevas etiquetas de los clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


### Guardar Resultados ###

#if (interactive()) {
#  saveRDS(pbmc, file = "pbmc3k_final.rds")
#}



### Información de la sesión de R ###
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
