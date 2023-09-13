########## Control de calidad ##########

## Paquetes de este capítulo
library("scRNAseq") # para descargar datos de ejemplo
library("AnnotationHub") # para obtener información de genes
library("scater") # para gráficas y control de calidad
library("BiocFileCache") # para descargar datos
library("DropletUtils") # para detectar droplets
library("Matrix") # para leer datos en formatos comprimidos


### Ejercicio: entendiendo addPerCellQC ###

## Obtención de los datos
library("scRNAseq")

sce.416b <- LunSpikeInData(which = "416b")

sce.416b$block <- factor(sce.416b$block)

## Descargar los archivos de anotación desde Ensembl vía AnnotationHub
library("AnnotationHub")

ah <- AnnotationHub()

query(ah, c("Mus musculus", "Ensembl", "v97")) # Búsqueda

## Obtener la posición del cromosoma para cada gen
ens.mm.v97 <- ah[["AH73905"]]

location <- mapIds(
  ens.mm.v97,
  keys = rownames(sce.416b),
  keytype = "GENEID",
  column = "SEQNAME"
)

## Identificar los genes mitocondriales
is.mito <- which(location == "MT")

library("scater")

## Añadir información adicional addPerCellQC en colData
sce.416b <- addPerCellQC(sce.416b,
                         subsets = list(Mito = is.mito)
)

## Generar una gráfica de boxplots del número de genes por bloque (block) de células
with(colData(sce.416b), boxplot(detected ~ block))


### Gráficas sobre medidas de control de calidad (QC) ###

## Grafico de dispersion de la relación entre los valores en "block" (x) y "detected" (y)
plotColData(sce.416b, x = "block", y = "detected")

## Grafico de dispersion de la relación entre los valores en "block" (x) y "detected" (y)
## aplicando una escala logarítmica
plotColData(sce.416b, x = "block", y = "detected") +
  scale_y_log10()

## Grafico de dispersion de la relación entre los valores en "block" (x) y "detected" (y)
## usando la variable "phenotype" para distinguir los puntos
## aplicando una escala logarítmica
## y dividido en paneles según el fenotipo ("phenotype")
plotColData(sce.416b,
            x = "block",
            y = "detected",
            other_fields = "phenotype"
) +
  scale_y_log10() +
  facet_wrap(~phenotype)
