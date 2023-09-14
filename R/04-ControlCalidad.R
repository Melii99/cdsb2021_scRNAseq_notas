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


### Ejercicio: gráficas QC ERCC ###

## Adaptar el código de las gráficas anteriores para otras variable de control de calidad

## Grafico de dispersion de la relación entre los valores en "block" (x) y "altexps_ERCC_percent" (y)
plotColData(sce.416b, x = "block", y = "altexps_ERCC_percent")

## Grafico de dispersion de la relación entre los valores en "block" (x) y "altexps_ERCC_percent" (y)
## usando la variable "phenotype" para distinguir los puntos
## y dividido en paneles según el fenotipo ("phenotype")
plotColData(sce.416b,
            x = "block",
            y = "altexps_ERCC_percent",
            other_fields = "phenotype"
) +
  facet_wrap(~phenotype)


### Eliminar células de baja calidad ###

## Valores de límite ejemplo
qc.lib <- sce.416b$sum < 100000
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10
qc.mito <- sce.416b$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

## Obtener un resumen del número de células eliminadas por cada filtro
DataFrame(
  LibSize = sum(qc.lib),
  NExprs = sum(qc.nexprs),
  SpikeProp = sum(qc.spike),
  MitoProp = sum(qc.mito),
  Total = sum(discard)
)

## Usando isOutlier() para determinar los valores de corte

qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower")

qc.nexprs2 <- isOutlier(sce.416b$detected,
                        log = TRUE,
                        type = "lower"
)

qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent,
                       type = "higher"
)

qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent,
                      type = "higher"
)

discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

## Extraer los límites de valores (thresholds)
attr(qc.lib2, "thresholds")

attr(qc.nexprs2, "thresholds")

## Obtener un resumen del número de células eliminadas por cada filtro
DataFrame(
  LibSize = sum(qc.lib2),
  NExprs = sum(qc.nexprs2),
  SpikeProp = sum(qc.spike2),
  MitoProp = sum(qc.mito2),
  Total = sum(discard2)
)

## Más pruebas
plotColData(sce.416b,
            x = "block",
            y = "detected",
            other_fields = "phenotype"
) +
  scale_y_log10() +
  facet_wrap(~phenotype)

## Determinar el bloque (batch) de muestras (fonotipo-bloque)
batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)

## Versión de isOutlier() que toma en cuenta los bloques de muestras
qc.lib3 <- isOutlier(sce.416b$sum,
                     log = TRUE,
                     type = "lower",
                     batch = batch
)

qc.nexprs3 <- isOutlier(sce.416b$detected,
                        log = TRUE,
                        type = "lower",
                        batch = batch
)

qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent,
                       type = "higher",
                       batch = batch
)

qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent,
                      type = "higher",
                      batch = batch
)

discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

## Extraer los límites de valores (thresholds)
attr(qc.lib3, "thresholds")

attr(qc.nexprs3, "thresholds")

## Obtener un resumen del número de células eliminadas por cada filtro
DataFrame(
  LibSize = sum(qc.lib3),
  NExprs = sum(qc.nexprs3),
  SpikeProp = sum(qc.spike3),
  MitoProp = sum(qc.mito3),
  Total = sum(discard3)
)

### Ejercicio: filtrado de células ###

# R. Sí fue necesario qc.lib para crear discard
# ¿Cúal filtro fue más estricto? ¿discard o discard2? R. discard
# Al considerar el grupo de cada muestra (batch), ¿descartamos más células usando un valor de límite automático? R. si
