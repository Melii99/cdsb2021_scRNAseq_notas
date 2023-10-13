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
