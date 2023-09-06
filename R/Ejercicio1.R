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
