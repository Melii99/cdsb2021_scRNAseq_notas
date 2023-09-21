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


### Ejercicio: ERCC Grun et al ###


## Adaptar el código de sce.416b para los datos de Grun et al y reproducir la imagen
# Fíjate en que variables de colData() estamos graficando.
# ¿Existe la variable discard en colData()?
# ¿Qué variable tiene valores de D10, D17, D2, D3 y D7?

## Agregar la información de las células filtradas (crear discard en colData)
sce.grun$discard <- discard.ercc2

## Grafica de ERCC_detected vs ERCC_sum
## Coloreada segun si es un outlier o no
## y paneles separados para cada valor único de "donor" (D10, D17, D2, D3 y D7)
plotColData(sce.grun,
            x = "altexps_ERCC_detected",
            y = "altexps_ERCC_sum",
            colour_by = "discard",
            other_fields = "donor"
            ) +
  facet_grid(~ donor)


### Identificando droplets vacíos con datos de PBMC ###

## Opciones algorítmicas para detecar los droplets vacíos ##

## Descargar los datos
library("BiocFileCache")

bfc <- BiocFileCache()
raw.path <-
  bfcrpath(
    bfc,
    file.path(
      "http://cf.10xgenomics.com/samples",
      "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
    )
  )

untar(raw.path, exdir = file.path(tempdir(), "pbmc4k")) # Extraer

## Leer los datos en R
library("DropletUtils")
library("Matrix")
library("scater")

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)
bcrank <- barcodeRanks(counts(sce.pbmc)) # Cálculo de ranking a partir de counts

sce.pbmc <- addPerCellQC(sce.pbmc) # añadir metricas de calidad

## Obtener los valores únicos (no repetidos)
uniq <- !duplicated(bcrank$rank)

## Graficar los UMIs en log. base 10
## indicando los valores de corte (knee y punto de inflexión)
plot(
  bcrank$rank[uniq],
  bcrank$total[uniq],
  log = "xy",
  xlab = "Rank",
  ylab = "Total UMI count",
  cex.lab = 1.2
)
abline(
  h = metadata(bcrank)$inflection,
  col = "darkgreen",
  lty = 2
)
abline(
  h = metadata(bcrank)$knee,
  col = "dodgerblue",
  lty = 2
)
legend(
  "bottomleft",
  legend = c("Inflection", "Knee"),
  col = c("darkgreen", "dodgerblue"),
  lty = 2,
  cex = 1.2
)

## Se pueden encontrar los droplets vacíos usando emptyDrops() !!!

## Usando DropletUtils para encontrar los droplets
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc)) # Detectar empty droplets y guardar objeto

#e.out (objeto emptyDrops) incluye FDR (false discovery rate)

## Revisar emptyDrops para una explicación de porque hay valores NA
summary(e.out$FDR <= 0.001) #Cuales tienen un FDR menor al .001


## R. La función asume que cuentas UMI menores al valor "lower" ya son droplets vacios
## y no calcula nada


## Otra vez usamos emptyDrops ahora especificando un límite menor
set.seed(100)
limit <- 100 # Especificar limite (menor)
# test.ambient=TRUE indica que queremos resultados también para los barcodes con valor igual o menor a "lower"
all.out <- emptyDrops(counts(sce.pbmc), lower = limit, test.ambient = TRUE)

# Idealmente, este histograma debería verse uniforme.
# Picos grandes cerca de cero indican que los _barcodes_
# con un número total de cuentas menor a "lower" no son
# de origen ambiental.
hist(all.out$PValue[all.out$Total <= limit &
                      all.out$Total > 0],
     xlab = "P-value",
     main = "",
     col = "grey80"
)

## ¿Por qué emptyDrops() regresa valores NA?
#R. La función no calcula nada para UMIs con valor igual o menor a "lower"
## ¿Los valores p son iguales entre e.out y all.out? R. No (NAs)
## ¿Son iguales si obtienes el subconjunto de valores que no son NA? R. Si


### Filtrado de expresión mitocondrial adicional ###

''' Falta obtener sce.pbmc$subsets_MT_percent

## Después de filtar los droplets, el filtrado por expresión mitocondrial ayuda
## a eliminar células de baja calidad.


## FDR menor o igual a 0.001 en e.out (al menos una célula)
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

## Identificar genes mitocondriales en sce.pbmc y almacena los índices  en is.mito
is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)

## Busqueda de outlayers (con isOutlier) para el subconjunto mitocondrial de sce.pbmc
discard.mito <- isOutlier(sce.pbmc$subsets_MT_percent, type = "higher")

## Gráficar del número total de lecturas y el porcentaje mitocondrial
plot(
  sce.pbmc$sum,
  sce.pbmc$subsets_MT_percent,
  log = "x",
  xlab = "Total count",
  ylab = "Mitochondrial %"
)
abline(h = attr(discard.mito, "thresholds")["higher"], col = "red")


### Ejercicio avanzado ###

# Volvamos a crear sce.pbmc para poder usar plotColData() y visualizar la relación entre total
# y los niveles de expresión mitocondrial (en porcentaje) separando lo que pensamos que son
# droplets vacíos y las células de acuerdo a los resultados que ya calculamos de emptyDrops().
# El resultado final se verá como en la siguiente imagen.

## Leer los datos en R
library("DropletUtils")
library("Matrix")
library("scater")

fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)
bcrank <- barcodeRanks(counts(sce.pbmc)) # Cálculo de ranking a partir de counts

## Obtener los valores únicos (no repetidos)
uniq <- !duplicated(bcrank$rank)

## Usando DropletUtils para encontrar los droplets
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc)) # Detectar empty droplets y guardar objeto

## Filtrar los NAs
e.out <- e.out[!is.na(e.out$FDR), ]

## Añadir metricas de calidad
sce.pbmc <- addPerCellQC(sce.pbmc)

# Grafica
plotColData(
  sce.pbmc,
  x = "total",
  y = "subsets_MT_percent",
  colour_by = "discard",
  other_fields = "phenotype"
) +
  facet_grid(~ sce.pbmc$is_cell)

'''

### Discusión ¿Conviene eliminar datos? ###

# Eliminemos las células de calidad baja
# al quedarnos con las columnas del objeto sce que NO
# queremos descartar (eso hace el !)
filtered <- sce.416b[, !discard2]
# Alternativamente, podemos marcar
# las células de baja calidad
marked <- sce.416b
marked$discard <- discard2

## ¿Cúal de estos objetos es más grande?
# R. marked es más grande, se le está añadiendo una nueva columna con discard2
# mientras filtered es aún más pequeño que sce.416b
## ¿Cúal prefieres usar? R. Es bueno usar el objeto marked para no perder información
# (en caso de tener suficiente memoria, claro)

### ExperimentSubset ###
# Permite hacer varios subconjuntos de datos y mantenerlos en un sólo objeto para trabajar
# con estos subsets (ej. original, filtrado, filtrado2... dentro de un mismo objeto)


### Explorando datos de forma interactiva con iSEE ###

## CRear un objeto sencillo de tipo RangedSummarizedExperiment
library("SummarizedExperiment")

## Crear los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6

## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)

## Juntar toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Explorar el objeto resultante
rse

## Explorar el objeto rse de forma interactiva
library("iSEE")

if (interactive()) {
  iSEE::iSEE(rse)
}

## Explorar el objeto sce.416b de forma interactiva
if (interactive()) {
  iSEE::iSEE(sce.416b, appTitle = "sce.416b")
}


### Información de la sesión de R ###
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
