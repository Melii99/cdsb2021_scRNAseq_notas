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

dense.size / sparse.size
