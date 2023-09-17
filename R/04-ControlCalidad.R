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
qc.lib <- sce.416b$sum < 100000 # menos de 100k lecturas
qc.nexprs <- sce.416b$detected < 5000 # menos de 5k genes
qc.spike <- sce.416b$altexps_ERCC_percent > 10 # más de 10% de ERCC
qc.mito <- sce.416b$subsets_Mito_percent > 10 # más de 10% lecturas mitocondriales
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

## Obtener un resumen del número de células eliminadas por cada filtro
DataFrame(
  LibSize = sum(qc.lib),
  NExprs = sum(qc.nexprs),
  SpikeProp = sum(qc.spike),
  MitoProp = sum(qc.mito),
  Total = sum(discard)
)

## Usando isOutlier() para determinar los valores de corte (scatter::isOutlier)

qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower") # Suma  de recuentos de expresión génica (lecturas)

qc.nexprs2 <- isOutlier(sce.416b$detected, # Genes detectados
                        log = TRUE,
                        type = "lower"
)

qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent, # Porcentaje de ERCC
                       type = "higher"
)

qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent, # Porcentaje de lecturas mitocondriales
                      type = "higher"
)

discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2 # Filtrado

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

qc.lib3 <- isOutlier(sce.416b$sum, # Suma  de recuentos de expresión génica
                     log = TRUE,
                     type = "lower",
                     batch = batch # Se agrupan los valores antes de calcular los valores atípicos
)

qc.nexprs3 <- isOutlier(sce.416b$detected, # Genes detectados
                        log = TRUE,
                        type = "lower",
                        batch = batch # Se agrupan los valores antes de calcular los valores atípicos
)

qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent, # Porcentaje de ERCC
                       type = "higher",
                       batch = batch # Se agrupan los valores antes de calcular los valores atípicos
)

qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent, # Porcentaje de lecturas mitocondriales
                      type = "higher",
                      batch = batch # Se agrupan los valores antes de calcular los valores atípicos
)

discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3 # Filtrado

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


### Datos de Grun et al ###

## Cargar set de datos
sce.grun <- GrunPancreasData()

## Añadir métricas de calidad
sce.grun <- addPerCellQC(sce.grun)

## Gráfico de dispersión
plotColData(sce.grun, x = "donor", y = "altexps_ERCC_percent")

## ¿Qué patrón revela esta gráfica?
# Todas las muestras de los donadores D17, D2 y D7 tienen un porcentaje similar de ERCC
# con outliers bien maracdos mientras los donadores D10 y D3 tienen una distribución de ERCC más continua
# Considerando que el contenido de ERCC debe ser bajo (porcentajes bajos), puede
# que haya problemas con las muestras provenientes de D10 y D3.

## ¿Cúal de las siguientes gráficas identifica mejor las células de baja calidad?

## isOutlier() puede ayudarnos cuando un grupo de muestras tuvo más problemas que el resto
discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
                          type = "higher",
                          batch = sce.grun$donor # Los outliers se evalúan dentro de cada grupo de "donor" por separado
)

discard.ercc2 <- isOutlier(
  sce.grun$altexps_ERCC_percent,
  type = "higher",
  batch = sce.grun$donor,
  subset = sce.grun$donor %in% c("D17", "D2", "D7") #  Solo se evalúan los donantes "D17", "D2", "D7"
)

## isOutlier() tomando en cuenta el batch (discard.ercc toma en cuenta todos los donadores)
plotColData(
  sce.grun,
  x = "donor",
  y = "altexps_ERCC_percent",
  colour_by = data.frame(discard = discard.ercc) # Colorear según si es outlier (considerando todos los batch)
)

## isOutlier() tomando en cuenta batch y muestras que fallaron (discard.ercc2 toma en cuenta solo donadores no problematicos)
plotColData(
  sce.grun,
  x = "donor",
  y = "altexps_ERCC_percent",
  colour_by = data.frame(discard = discard.ercc2) # Colorear según si es outlier (considerando solo D17, D2, D7"
)

# La segunda opción, tomando en cuenta unicamente las muestras "no problemáticas"
# (solo de los donadores C17, D2 y D7) para isOutlayer() permite calcular mejor los
# verdaderos outlayers, pues el método asume que la mayoría de las células son de
# buena calidad, y al usar también a los donadores D10 y D3 obtenemos un corte más laxo.


### Gráficas de QC extra ###

## Otras grafica que se pueden hacer

## Agregar información sobre que células que tienen valores extremos (outliers)
sce.416b$discard <- discard2

## Hacer una gráfica para cada medida de control de calidad (QC)

## Gráfica QC de las lecturas detectadas de sce.416b por bloque y por fenotipo
## coloreada según si son o no valores extremos
## aplicando una escala logarítmica
plotColData(
  sce.416b,
  x = "block",
  y = "sum",
  colour_by = "discard",
  other_fields = "phenotype"
) +
  facet_wrap(~phenotype) +
  scale_y_log10()

## Gráfica QC de los genes detectados de sce.416b por bloque y por fenotipo
## coloreada según si son o no valores extremos
## aplicando una escala logarítmica
plotColData(
  sce.416b,
  x = "block",
  y = "detected",
  colour_by = "discard",
  other_fields = "phenotype"
) +
  facet_wrap(~phenotype) +
  scale_y_log10()

## Gráfica QC del porcentaje de ERCC de sce.416b por bloque y por fenotipo
## coloreada según si son o no valores extremos
plotColData(
  sce.416b,
  x = "block",
  y = "altexps_ERCC_percent",
  colour_by = "discard",
  other_fields = "phenotype"
) +
  facet_wrap(~phenotype)

## Gráfica QC del porcentaje de lecturas mitocondriales de sce.416b por bloque y por fenotipo
## coloreada según si son o no valores extremos
plotColData(
  sce.416b,
  x = "block",
  y = "subsets_Mito_percent",
  colour_by = "discard",
  other_fields = "phenotype"
) +
  facet_wrap(~phenotype)


## Otra gráfica de diagnóstico útil
plotColData(
  sce.416b,
  x = "sum",
  y = "subsets_Mito_percent",
  colour_by = "discard",
  other_fields = c("block", "phenotype")
) +
  facet_grid(block ~ phenotype)

















