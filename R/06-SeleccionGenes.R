########## Selección de Genes ##########

## Usualmente usamos datos scRNA-seq para caracterizar la heterogeneidad entre células
## Para hacer esto, usamos métodos como el clustering y la reducción de dimensionalidad
## Esto involucra resumir las diferencias por gen en una sola medida de (dis)similitud
## entre un par de células
## ¿Cuáles genes deberíamos usar para calcular esta medida de (dis)similitud?

## Deseamos seleccionar los genes altamente variables (High Variable Genes HVGs).
## Genes con una variación incrementada en comparación con otros genes que están
##siendo afectados por ruido técnico u otra variación biológica que no es de interés.


### Descarga de datos ###

library(BiocFileCache)

## Gestionar el almacenamiento en caché de archivos biológicos
bfc <- BiocFileCache()

## Descarga del archivo
raw.path <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples",
  "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
))

## Descomprimir el archivo
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
library(Matrix)

## Ruta al directorio donde se encuentran los datos extraídos
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")

## Cargar los datos de expresión desde fname en un objeto SingleCellExperiment (sce.pbmc)
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

sce.pbmc


###  Anotación ###

library(scater)

## Nombres únicos de filas (rownames) en función de las columnas "ID" y "Symbol"
## de los metadatos de fila (rowData)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol
)

##  base de datos de anotación genómica de Homo sapiens
library(EnsDb.Hsapiens.v86)

## Mapeo genómico
location <- mapIds(EnsDb.Hsapiens.v86,
                   keys = rowData(sce.pbmc)$ID,
                   column = "SEQNAME", keytype = "GENEID"
)



### Control de calidad QC ###

## Eliminar _droplets_ vacíos
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc)) # Detectar roplets vacíos
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)] # Retener solo las células que pasaron el filtro

## Calcular métricas de control de calidad (QC)
stats <- perCellQCMetrics(sce.pbmc,
                          subsets = list(Mito = which(location == "MT"))
)
## Detectar alta expresión de genes mitocondriales (y guardar índices)
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type = "higher"
)
## Eliminar las células con alta expresión de genes mitocondriales
sce.pbmc <- sce.pbmc[, !high.mito]


### Normalización de los datos ###
library(scran)
set.seed(1000)

## Normalización por decircunvolución (deconvolution)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)

## Transformación logarítmica
sce.pbmc <- logNormCounts(sce.pbmc)


### Preguntas de repaso ###

## ¿Cómo determinamos cuáles eran los genes mitocondriales?
## R. A partir de la anotación usando EnsDb.Hsapiens.v86, y después seleccionando
## los genes correspondientes (location == "MT")
## ¿Cómo decidimos filtrar las células?
## Para los resultados de emptyDrops() se fijó un límite de 0.1% FDR
## y para el filtro de genes mito. se usó isOutlier-"higher", el cual
## utiliza como límite 3 desviaciones sobre la mediana (MAD).
## ¿Puedes explicar como normalizamos los datos?
## Generando clusters (con quickCluster) y usandolos para calcular los factores
## de normalización correspondientes (computeSumFactors)
