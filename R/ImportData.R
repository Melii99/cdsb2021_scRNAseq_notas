########## Import Data ##########


### CellRanger ###

## Descargar los datos de ejemplo procesados con CellRanger

# Nota: A usar BiocFileCache solo tenemos que descargar los datos una vez.
library("BiocFileCache")

bfc <- BiocFileCache()

## Crear una URL de descarga para el archivo tar.gz
pbmc.url <-
  paste0(
    "http://cf.10xgenomics.com/samples/cell-vdj/",
    "3.1.0/vdj_v1_hs_pbmc3/",
    "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
  )

##  bfcrpath para descargar el archivo ZIP desde la URL
pbmc.data <- bfcrpath(bfc, pbmc.url)

## Extraer los archivos en un directorio temporal
untar(pbmc.data, exdir = tempdir())

## Enumerar los archivos que descargamos y que extrajimos
pbmc.dir <- file.path(
  tempdir(),
  "filtered_feature_bc_matrix" # Estos son los archivos típicos de CellRanger
)

list.files(pbmc.dir)

## Importar los datos como un objeto de tipo SingleCellExperiment
library("DropletUtils")

sce.pbmc <- read10xCounts(pbmc.dir)

## Revisar el objeto que acabamos de construir
sce.pbmc

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce.pbmc) / 1024^2

## Almacenar la información de CITE-seq como un experimento alternativo
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)

## Revisar el objeto sce.pbmc actualizado
sce.pbmc

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce.pbmc) / 1024^2


### scPipe ###

## Descargar datos de ejemplo procesados con scPipe
library("BiocFileCache")

bfc <- BiocFileCache()

sis_seq.url <-
  "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"

##  bfcrpath para descargar el archivo ZIP desde la URL
sis_seq.data <- bfcrpath(bfc, sis_seq.url)

## Extraer los archivos en un directorio temporal
unzip(sis_seq.data, exdir = tempdir())

## Enumerar algunos de los archivos que descargamos y extrajimos
sis_seq.dir <- file.path(
  tempdir(),
  "SIS-seq_script-master", # Estos son los archivos típicos de scPipe
  "data",
  "BcorKO_scRNAseq",
  "RPI10"
)

list.files(sis_seq.dir)

## Importar los datos como un objeto de tipo SingleCellExperiment
library("scPipe")

sce.sis_seq <- create_sce_by_dir(sis_seq.dir)

## Revisar el objeto que acabamos de construir
sce.sis_seq

## ¿Qué tan grande es el objeto de R (En MB)?
lobstr::obj_size(sce.sis_seq) / 1024^2


### Archivos Varios ###


#La dirección https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5522/E-MTAB-5522.processed.1.zip
#ya no se encuentra disponible!

'''

## Ejemplo de descarga de archivos varios
library("BiocFileCache")

bfc <- BiocFileCache()

## Archivo de cuentas
## Ya no se encuentra disponible :c

lun_counts.url <-
  paste0(
    "https://www.ebi.ac.uk/arrayexpress/files/",
    "E-MTAB-5522/E-MTAB-5522.processed.1.zip"
  )

##  bfcrpath para descargar el archivo ZIP desde la URL
lun_counts.data <- bfcrpath(bfc, lun_counts.url)

## Archivo colData
lun_coldata.url <-
  paste0(
    "https://www.ebi.ac.uk/arrayexpress/files/",
    "E-MTAB-5522/E-MTAB-5522.sdrf.txt"
  )

##  bfcrpath para descargar el archivo ZIP desde la URL
lun_coldata.data <- bfcrpath(bfc, lun_coldata.url)

## Extraer los archivos en un directorio temporal
lun_counts.dir <- tempfile("lun_counts.")
unzip(lun_counts.data, exdir = lun_counts.dir)

## Enumerar los archivos que descargamos y extrajimos
list.files(lun_counts.dir)


## Leer la matriz de cuentas (para una placa)
lun.counts <- read.delim(
  file.path(lun_counts.dir, "counts_Calero_20160113.tsv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

## Almacenar la información de la longitud de los genes para después
gene.lengths <- lun.counts$Length

## Convertir los datos de cuentas de genez a una matriz (quitamos las longitudes)
lun.counts <- as.matrix(lun.counts[, -1])


## Leer la información de las muestras (células)
lun.coldata <- read.delim(lun_coldata.data,
                          check.names = FALSE,
                          stringsAsFactors = FALSE
)

library("S4Vectors")

lun.coldata <- as(lun.coldata, "DataFrame")

## Poner en orden la información de las muestras para que sea idéntico al orden en la matriz de cuentas
m <- match(
  colnames(lun.counts),
  lun.coldata$`Source Name`
)

lun.coldata <- lun.coldata[m, ]

## Construir la tabla de información de los genes
lun.rowdata <- DataFrame(Length = gene.lengths)

## Construir el objeto de SingleCellExperiment
lun.sce <- SingleCellExperiment(
  assays = list(assays = lun.counts),
  colData = lun.coldata,
  rowData = lun.rowdata
)

## Revisar el objeto que acabamos de construir
lun.sce

'''

### Información de la sesión de R ###
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
