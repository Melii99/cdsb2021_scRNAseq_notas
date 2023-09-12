########## Import Data ##########


### CellRanger ###

## Descargar los datos de ejemplo procesados con CellRanger

# Nota: A usar BiocFileCache solo tenemos que descargar los datos una vez.
library("BiocFileCache")

bfc <- BiocFileCache()

## Crear una URL de descarga para el archivo tar.gz
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
