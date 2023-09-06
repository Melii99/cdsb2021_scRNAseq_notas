########## Ejercicio 1 ##########

### Data ###

## Obtener el set de datos de 416b
library("scRNAseq")
sce.416b <- LunSpikeInData(which = "416b")

## Carga el paquete SingleCellExperiment
library("SingleCellExperiment")

#r# Extrer la matriz de cuentas del set de datos de 416b
counts.416b <- counts(sce.416b)

## Clase y dimensiones de la matriz de cuentas
class(counts.416b)
dim(counts.416b)


### SCE object ###

## Construir un nuevo objeto SCE de la matriz de cuentas
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

## Revisa el objeto que acabamos de crear
sce

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2

## Accesar la matriz de cuenta del compartimento (slot) "assays"
# assays(sce, "counts")
# OJO: ¡esto puede inundar tu sesión de R!

# 1. El método general
#assay(sce, "counts")[110:115, 1:3] # gene, cell

# 2. El método específico para accesar la matriz de cuentas "counts"
#counts(sce)[110:115, 1:3]


### Agregar más assays ###

## Usar logNormCounts para normalizar
sce <- scater::logNormCounts(sce)

## Revisa el objeto actualizado (assay adicional logcounts)
sce

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2

## Método general para accesar la matriz de cuentas "logcounts"
#assay(sce, "logcounts")[110:115, 1:3]

## Método específico para accesar la matriz de cuentas transformadas "logcounts"
#logcounts(sce)[110:115, 1:3]

### Agregar un assay adicional de manera manual ###

## Suma 100 a counts y se agrega como un nuevo assay "counts_100"
assay(sce, "counts_100") <- assay(sce, "counts") + 100

## Enumerar los "assays" de sce (numero y nombre de los assays)
assays(sce)

## Nombres de los assays
assayNames(sce)


### Metadata (colData) ###

## Extraer la información de las muestras (metadata) del set de datos de 416b
colData.416b <- colData(sce.416b)

## Explorar datos
table(colData.416b$phenotype)
## fue en varios dias?
table(colData.416b$block)

## Agregar algo de esa información (fenotipo y block) a nuestro objeto SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]

## Revisar el objeto actualizado
sce

## Accesar a la información de las muestras (metadata) en nuestro SCE
colData(sce)

## Accesar una columna específica de la información de las muestras (metadata)
table(sce$block)

## Otra manera de accesar una columna específica de metadata
table(colData(sce)$block)

### Agregar columnas nuevas a colData: Ej. addPerCellQC ###

## Añadir datos de control de calidad
sce <- scater::addPerCellQC(sce.416b)

## Accesar a la metadata (información de las muestras) del objeto SCE actualizado
colData(sce)

## Revisar el objeto SCE actualizado
sce

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2

## Agrega las cuentas normalizadas (logNormCounts) de nuevo
sce <- scater::logNormCounts(sce)

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2

## Ejemplo: obtener el subconjunto de células de fenotipo "wild type"
## Recordatorio: las células son columnas del SCE
sce[, sce$phenotype == "wild type phenotype"]


### metadata de features (rowData) ###

## Accesar la información de los genes de  SCE (actualmente vacío)
rowData(sce)

### Agregar columnas nuevas a rowData: Ej. addPerFeatureQC ###

sce <- scater::addPerFeatureQC(sce)

## Accesar a la metadata (información de las muestras) del SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2


### Descargar los archivos de anotación  ###

## Archivos para anotación de la base de datos de Ensembl vía AnnotationHub
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

## Obtener la posición del cromosoma para cada gen
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
                     keys = rownames(sce),
                     keytype = "GENEID",
                     column = "SEQNAME"
)
rowData(sce)$chromosome <- chromosome

## Accesar a la metadata (información de las muestras) en el SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2

## Ejemplo: obtener el subconjunto de datos donde los genes están en el cromosoma 3
## NOTA: which() fue necesario para lidear con los nombres de cromosoma que son NA
sce[which(rowData(sce)$chromosome == "3"), ]


### Slot Metadata ###

## Accesar a la información de nuestro experimento usando metadata() (actualmente vacía)
metadata(sce)

## Es posible agregar cualquier tipo de información a metadata()
## Ej. Agregar genes favoritos y analista
metadata(sce) <- list(
  favourite_genes = c("Shh", "Nck1", "Diablo"),
  analyst = c("Pete")
)

## Accesar a la información de metadata() de SCE actualizado
metadata(sce)


### Slot de reducción de dimensiones ###

## Agregarlos componentes principales (PCs) de las logcounts
sce <- scater::runPCA(sce)

## Revisa el objeto actualizado (PCA agregado a reducedDimNames)
sce

## Accesar a la matriz de PCA del componente (slot) reducedDims (6 cells, 3 PCs)
reducedDim(sce, "PCA")[1:6, 1:3]

## Agregar una representación de los logcounts en t-SNE
sce <- scater::runTSNE(sce)

# Revisa el objeto actualizado (TSNE agregado a reducedDimNames)
sce

## Accesar a la matriz de t-SNE en el componente (slot) de reducedDims
head(reducedDim(sce, "TSNE"))


## Agregar manualmente una representación de los logcounts en UMAP
u <- uwot::umap(t(logcounts(sce)), n_components = 2)
## Agrega la matriz de UMAP al componente (slot) reducedDims
reducedDim(sce, "UMAP") <- u
# Accesa a la matriz de UMAP desde el componente (slot) reducedDims
head(reducedDim(sce, "UMAP"))

## Enumerar los resultados de reducción de dimensiones de SCE
reducedDims(sce)


### Alternative experiments ###

## Extraer la información de ERCC de nuestro SCE para el set de datos de 416b
ercc.sce.416b <- altExp(sce.416b, "ERCC")

## Inspeccionar el SCE para los datos de ERCC
ercc.sce.416b

## Agregar el SCE de ERCC como un experimento alternativo a nuestro SCE
altExp(sce, "ERCC") <- ercc.sce.416b

## Revisar nuestro objeto SCE actualizado (altExpNames - ERCC SIRV)
sce

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2

## Enumerar los experimentos alternativos almacenados en nuestro objeto sce
altExps(sce)


## Si creamos un subconjunto del SCE por muestra (célula)
sce.subset <- sce[, 1:10]
ncol(sce.subset)
## automáticamente se obtiene el subconjunto de los experimentos alternativos
ncol(altExp(sce.subset))

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce) / 1024^2


### Size factors ###

## Extraer los factores de tamaño (size factors)
# Estos fueron añadidos a nuestro objeto cuando corrimos scater::logNormCounts(sce)
head(sizeFactors(sce))

## "Automáticamente" reemplaza los factores de tamaño
sce <- scran::computeSumFactors(sce)
head(sizeFactors(sce))

## "Manualmente" reemplaza los factores de tamaño
sizeFactors(sce) <- scater::librarySizeFactors(sce)
head(sizeFactors(sce))
