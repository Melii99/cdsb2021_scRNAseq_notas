########## Import Data ##########

### CellRanger ###

## Descargar los datos de ejemplo procesados con CellRanger

# Nota: A usar BiocFileCache solo tenemos que descargar los datos una vez.
library("BiocFileCache")

bfc <- BiocFileCache()
pbmc.url <-
  paste0(
    "http://cf.10xgenomics.com/samples/cell-vdj/",
    "3.1.0/vdj_v1_hs_pbmc3/",
    "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
  )

pbmc.data <- bfcrpath(bfc, pbmc.url)

## Extraer los archivos en un directorio temporal
untar(pbmc.data, exdir = tempdir())

## Enumerar los archivos que descargamos y que extrajimos
pbmc.dir <- file.path(
  tempdir(),
  "filtered_feature_bc_matrix"
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
