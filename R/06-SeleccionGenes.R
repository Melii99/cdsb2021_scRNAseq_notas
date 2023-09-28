########## Selección de Genes ##########

## Usualmente usamos datos scRNA-seq para caracterizar la heterogeneidad entre células
## Para hacer esto, usamos métodos como el clustering y la reducción de dimensionalidad
## Esto involucra resumir las diferencias por gen en una sola medida de (dis)similitud
## entre un par de células
## ¿Cuáles genes deberíamos usar para calcular esta medida de (dis)similitud?

## Deseamos seleccionar los genes altamente variables (High Variable Genes HVGs).
## Genes con una variación incrementada en comparación con otros genes que están
##siendo afectados por ruido técnico u otra variación biológica que no es de interés.


### Descarga de datos ###

library(BiocFileCache)

## Gestionar el almacenamiento en caché de archivos biológicos
bfc <- BiocFileCache()

## Descarga del archivo
raw.path <- bfcrpath(bfc, file.path(
  "http://cf.10xgenomics.com/samples",
  "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
))

## Descomprimir el archivo
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
library(Matrix)

## Ruta al directorio donde se encuentran los datos extraídos
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")

## Cargar los datos de expresión desde fname en un objeto SingleCellExperiment (sce.pbmc)
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

sce.pbmc


###  Anotación ###

library(scater)

## Nombres únicos de filas (rownames) en función de las columnas "ID" y "Symbol"
## de los metadatos de fila (rowData)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol
)

##  base de datos de anotación genómica de Homo sapiens
library(EnsDb.Hsapiens.v86)

## Mapeo genómico
location <- mapIds(EnsDb.Hsapiens.v86,
                   keys = rowData(sce.pbmc)$ID,
                   column = "SEQNAME", keytype = "GENEID"
)



### Control de calidad QC ###

## Eliminar _droplets_ vacíos
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc)) # Detectar roplets vacíos
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)] # Retener solo las células que pasaron el filtro

## Calcular métricas de control de calidad (QC)
stats <- perCellQCMetrics(sce.pbmc,
                          subsets = list(Mito = which(location == "MT"))
)
## Detectar alta expresión de genes mitocondriales (y guardar índices)
high.mito <- isOutlier(stats$subsets_Mito_percent,
                       type = "higher"
)
## Eliminar las células con alta expresión de genes mitocondriales
sce.pbmc <- sce.pbmc[, !high.mito]


### Normalización de los datos ###
library(scran)
set.seed(1000)

## Normalización por decircunvolución (deconvolution)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)

## Transformación logarítmica
sce.pbmc <- logNormCounts(sce.pbmc)


### Preguntas de repaso ###

## ¿Cómo determinamos cuáles eran los genes mitocondriales?
## R. A partir de la anotación usando EnsDb.Hsapiens.v86, y después seleccionando
## los genes correspondientes (location == "MT")

## ¿Cómo decidimos filtrar las células?
## R. Para los resultados de emptyDrops() se fijó un límite de 0.1% FDR
## y para el filtro de genes mito. se usó isOutlier-"higher", el cual
## utiliza como límite 3 desviaciones sobre la mediana (MAD).

## ¿Puedes explicar como normalizamos los datos?
## R. Generando clusters (con quickCluster) y usandolos para calcular los factores
## de normalización correspondientes (computeSumFactors)


### Dataset ilustrativo: 416B ###

## Línea celular de células mieloides progenitoras inmortalizadas de ratón usando SmartSeq2

## Descarga de datos
library(scRNAseq)

sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

## Anotación
library(AnnotationHub)

ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97,
                                   keys = rownames(sce.416b),
                                   keytype = "GENEID", column = "SYMBOL"
)

rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97,
                                    keys = rownames(sce.416b),
                                    keytype = "GENEID", column = "SEQNAME"
)

library(scater)

rownames(sce.416b) <- uniquifyFeatureNames(
  rowData(sce.416b)$ENSEMBL,
  rowData(sce.416b)$SYMBOL
)

## Control de calidad (QC)
mito <- which(rowData(sce.416b)$SEQNAME == "MT")

stats <- perCellQCMetrics(sce.416b, subsets = list(Mt = mito))

qc <- quickPerCellQC(stats,
                     percent_subsets = c("subsets_Mt_percent", "altexps_ERCC_percent"),
                     batch = sce.416b$block
)

sce.416b <- sce.416b[, !qc$discard]

## Normalización
library(scran)

sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)

### Preguntas de repaso ###

## ¿Cómo determinamos cuáles eran los genes mitocondriales?
## R. A partir de la anotación usando ens.mm.v9, y después seleccionando
## los genes correspondientes (SEQNAME == "MT")

## ¿Cómo decidimos filtrar las células?
## R. Se generaron  QC metrics (perCellQCMetrics) para el subset de Mito. y después
## obtener los porcentajes de expresión génica de mito y ERCC (quickPerCellQC) para
## finalmente filtrar a partir de esos resultados

## ¿Puedes explicar como normalizamos los datos?
## R. Por escalamiento, utilizando computeSumFactors para obtener los factores de
## normalización y haciendo una transformación logaritmica con logNormCounts



### Cuantificando la varianza por gen ###

## Varianza de los log-counts ##

## Selección del feature basado en los log-counts (usadas en análisis posteriores)
## La transformación log no logra la estabilización de la varianza perfecta,
## así que se requiere modelar la relación de la varianza-media de los features


## Enfoque simple: Calcular la varianza de los log-counts para cada gen (ignorando
## grupos experimentales) y ordenar los genes del más-al-menos variable

## Un enfoque más sofisticado: Se calcula la varianza de los log-counts para cada gen
## (ignorando grupos experimentales), se modela la relación entre la media y la varianza
## de los log-counts para estimar la variación técnica. Se estima la varianza
## biológica sustrayendo la varianza técnica de la varianza total y se ordenan
## los genes de la variable de mayor-a-menor biológicamente

## Cálculo de la varianza con scran::modelGeneVar ##

## Supuestos:
## La abundancia de los perfiles de expresión de la mayoría de los genes están
## dominados por el ruido aleatorio técnico.
## Una tendencia representa un estimado del ruido técnico como una función de la abundancia
## Así, es posible descomponer la varianza total de cada gen en un componente técnico y uno biológico
## Genes con una gran varianza biológica son considerados interesantes

## Varianza de las log-counts con scran::modelGeneVar
library(scran)
dec.pbmc <- modelGeneVar(sce.pbmc)

## Visualizaciíon de la relación entre la media y la varianza
fit.pbmc <- metadata(dec.pbmc)

plot(fit.pbmc$mean, fit.pbmc$var,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression"
)

curve(fit.pbmc$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


### Ejercicios ###

## ¿Qué tipo de objeto nos regresó modelGeneVar()? R. DFrame

## ¿dec.pbmc es una tabla? ¿O contiene mayor información? R. Es un data frame con
## 33694 filas y 6 columnas y contiene más información dentro de metadata(dec.pbmc)

## ¿Qué tipo de objeto es fit.pbmc y que objetos con nombres contiene?
## R. Es una lista, y contiene "mean", "var", "trend" y "std.dev"

## ¿Qué tipo de objeto es fit.pbmc$trend? R. función

## ¿Donde podemos encontrar más detalles de esta función? R. ?fitTrendVar


### Ordenando por genes más interesantes ###

## Ordenamos los genes por mayor "varianza biológica" (más interensantes)
dec.pbmc[order(dec.pbmc$bio, decreasing = TRUE), ]



### Coeficiente de variación de las cuentas ###

## El coeficiente de variación de las cuentas al cuadrado (CV2) es una alternativa
## a la varianza de los log-counts (uso de desviación estándar)

## Se calcula usando las cuentas en lugar de los log-counts !!!
## CV es el ratio de la desviación estándar a la media y está muy relacionada con
## el parámetro de dispersión de la distribución binomial negativa (edgeR y DESeq2)

## Cálculo del Coeficiente de variación con modelGeneCV2 ##

## Modela la relación de la media de la varianza cuando se considera la relevancia de cada gen
## Asume que la mayoría de los genes contienen ruido aleatorio y que la tendencia
## captura la mayoría de la variación técnica
## Genes con un gran CV2 que se desvían fuertemente de la tendencia es probable que
## representen genes afectados por la estructura biológica
## Usa la proporción (en lugar de la diferencia) del CV2 a la tendencia


## Coeficiente de variación con modelGeneCV2
dec.cv2.pbmc <- modelGeneCV2(sce.pbmc)

## Visuualizar la relación con la media
fit.cv2.pbmc <- metadata(dec.cv2.pbmc)

plot(fit.cv2.pbmc$mean, fit.cv2.pbmc$cv2,
     log = "xy"
)

curve(fit.cv2.pbmc$trend(x),
      col = "dodgerblue",
      add = TRUE, lwd = 2
)

### Ordenando los genes por coeficiente de variación ###

## Ordenamos los genes "más interesantes" por su coeficiente de variación
dec.cv2.pbmc[order(dec.cv2.pbmc$ratio, decreasing = TRUE), ]


### Varianza de los log-counts vs coeficiente de variación ###

## Generalmente se usa la varianza de los log-counts !!!

## ## CV2 tiende a tener otorgar rangos altos a genes altamente variables (HGVs)
## con bajos niveles de expresión
## Éstos son dirigidos por una sobreregulación en subpoblaciones raras
## Puede asignar un alto rango a genes que no son de nuestro interés con varianza baja absoluta
## La variación descrita por el CV2 de las cuentas es menos relevante para los
## procedimientos que operan en los log-counts



### Cuantificando el ruido técnico ###

## Aunque ajustamos una línea de tendencia a todos los genes endógenos y asumimos
## que la mayoría de los genes no están dominados por ruido técnico, en la práctica,
## todos los genes expresados tienen algún nivel de variabilidad biológica

## Es mejor que pensemos estos estimados como una variación técnica más la
## variación biológica que no es interesante

## Podemos revisar dos enfoques:
## Cuando tenemos spike-ins
## Cuando no tenemos spike-ins


### Cuantificando el ruido técnico en presencia de spike-ins ###

## scran::modelGeneVarWithSpikes Ajusta la tendencia solo con los spike-ins
## (que deberían estar afectados solamente por variación técnica)

## Modelo de la varianza de la expresión de los genes en relación con los ERCC
dec.spike.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC")

## Ordenar los genes por "más interesantes" - mayor varianza
dec.spike.416b[order(dec.spike.416b$bio, decreasing = TRUE), ]

## Visualizar la varianza de la expresión de los genes en relación con los ERCC
plot(dec.spike.416b$mean, dec.spike.416b$total,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression"
)

fit.spike.416b <- metadata(dec.spike.416b)

points(fit.spike.416b$mean, fit.spike.416b$var,
       col = "red", pch = 16
)

curve(fit.spike.416b$trend(x),
      col = "dodgerblue",
      add = TRUE, lwd = 2
)


### Cuantificando el ruido técnico en ausencia de spike-ins ###

## scran::modelGeneVarByPoisson Realiza asunciones estadísticas acerca del ruido
## Las cuentas UMI típicamente muestran una variación cercana a Poisson si solo
## consideramos ruido técnico de la preparación de las librerías y la secuenciación

## modelGeneVarByPoisson() realiza simulaciones, por lo que necesitamos “ajustar
## la “semilla” para obtener resultados reproducibles

## modelGeneVarByPoisson() pueden también simular una variación binomial negativa
## (variación de Poisson sobredispersada)

## Cuantificando la varianza (cuantif. ruido técnico) asumiendo una distribución poisson
set.seed(0010101)
dec.pois.pbmc <- modelGeneVarByPoisson(sce.pbmc)

## Ordenar los genes por "más interesantes" - mayor varianza
dec.pois.pbmc[order(dec.pois.pbmc$bio, decreasing = TRUE), ]

## Visualizar la varianza de la expresión de los genes obtenido por modelGeneVarByPoisson
plot(dec.pois.pbmc$mean, dec.pois.pbmc$total,
     pch = 16, xlab = "Mean of log-expression",
     ylab = "Variance of log-expression"
)

curve(metadata(dec.pois.pbmc)$trend(x),
      col = "dodgerblue", add = TRUE
)



### Recordando propiedades de los datos de sce.416b ###

## línea celular de células inmortalizadas mieloides progenitoras de ratón con SmartSeq2

## cantidad constante de spike-in ERCC RNA se agregó a cada lisado celular antes
## de la preparación de la librería


### Considerando factores experimentales ###

## Los datos que contienen múltiples batches, muy seguido presentan efecto de bloque
## que pueden crear HGVs (genes altamente variables) artificiales !!!

## Se debe identificar los HGVs en cada batch y combinarlos en una única lista de HGVs

## Verificar que existen batches ("block")
names(colData(sce.416b))

## Modelo de la varianza de la expresión de los genes en relación con los ERCC
## con  modelGeneVarWithSpikes por bloque (evitar HVGs artificiales por batch effect)
dec.block.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC",
                                         block = sce.416b$block) # Usar parámetro block

## Ordenar los genes por "más interesantes" - mayor varianza
dec.block.416b[order(dec.block.416b$bio, decreasing = TRUE), ]




###  Seleccionando genes altamante variables (high-variable genes, HVGs) ###


## Hasta ahora hemos ordenado los genes del más al menos interesantemente variable
## ¿Qué tanto debemos de bajar en la lista para seleccionar nuestros HVGs?

## Es difícil determinar el balance óptimo porque el rudio en un contexto podría
## ser una señal útil en otro contexto, por ello existen varias estrategias:


## Seleccionando HVGs sobre la métrica de varianza (bio) ##

## La estrategia más simple es seleccionar los top-X genes con los valores más grandes
## para la métrica relevante de varianza (Ej. top values of scran::modelGeneVar())

## Pro: El usuario puede controlar directamente el número de HVGs
## Contra: ¿Qué valor de X se debe usar?

## Usar getTopHVGs para optener los genes más relevantes ##

## Obteniendo los top 1000 genes con la mayor varianza encontrados por modelGeneVar()
hvg.pbmc.var <- getTopHVGs(dec.pbmc, n = 1000) # Works with modelGeneVar() output
str(hvg.pbmc.var)

## Obteniendo los top 1000 genes con la mayor varianza encontrados por modelGeneVarWithSpikes()
hvg.416b.var <- getTopHVGs(dec.spike.416b, n = 1000) # Works with modelGeneVarWithSpikes() output
str(hvg.416b.var)

## Obteniendo los top 1000 genes con la mayor varianza encontrados por
## modelGeneVarWithSpikes() por batches
hvg.pbmc.cv2 <- getTopHVGs(dec.cv2.pbmc, # Also works with modelGeneCV2()
                           var.field = "ratio", n = 1000 # but note `var.field`
)
str(hvg.pbmc.cv2)

## Estrategias para seleccionar X

## Podemos asumir que, por ejemplo, 5% de los genes están diferencialmente expresados
## Establece X como el 5% de los genes

## Normalmente no conocemos el número de genes diferencialmente expresados desde antes,
## por lo tanto, solo hemos cambiado un número arbitrario por otro número arbitrario


### Seleccionando HVGs de acuerdo a su significancia estadística ###


## Establecer un límite fijo en alguna métrica de significancia estadística.
## Ej. p-value para cada gen: seleccionar todos los genes con un p-valor ajustado menor que 0.05

## Recordatorio: las pruebas estadísticas siempre dependen del tamaño de la muestra

## Pros: * Fácil de implementar * Menos predecible que la estrategia de los top-X
## Contras: * Podría priorizar genes con significancia estadística fuerte
## en vez de significancia biológica fuerte

## Obteniendo los genes con un p-value ajustado menor a 0.05 encontrados por modelGeneVar()
hvg.pbmc.var.2 <- getTopHVGs(dec.pbmc, fdr.threshold = 0.05) # Works with modelGeneVar() output
str(hvg.pbmc.var.2)

## Obteniendo los genes con un p-value ajustado menor a 0.05 encontrados por modelGeneVarWithSpikes()
hvg.416b.var.2 <- getTopHVGs(dec.spike.416b, # Works with modelGeneVarWithSpikes() output
                             fdr.threshold = 0.05
)
str(hvg.416b.var.2)

## Obteniendo los genes con un p-value ajustado menor a 0.05 encontrados por
##  modelGeneVarWithSpikes() por batches
hvg.pbmc.cv2.2 <- getTopHVGs(dec.cv2.pbmc, # Also works with modelGeneCV2()
                             var.field = "ratio", fdr.threshold = 0.05 # but note `var.field`
)
str(hvg.pbmc.cv2.2)



### Seleccionando genes por arriba de la tendencia media-varianza ###

## Selecciona todos los genes con una varianza biológica positiva

## Este es un extremo del equilibrio sesgo-varianza que minimiza el sesgo con el
## costo de maximizar el ruido

## Funciona mejor si tenemos datasets altamente heterogeneos que contienen
## muchos tipos celulares diferentes

## Obteniendo los genes por arriba de la tendencia media-varianza encontrados por modelGeneVar()
hvg.pbmc.var.3 <- getTopHVGs(dec.pbmc, var.threshold = 0) # Works with modelGeneVar() output
str(hvg.pbmc.var.3)

## Obteniendo los genes por arriba de la tendencia media-varianza encontrados por modelGeneVarWithSpikes()
hvg.416b.var.3 <- getTopHVGs(dec.spike.416b, # Works with modelGeneVarWithSpikes() output
                             var.threshold = 0
)
str(hvg.416b.var.3)

## Obteniendo los genes por arriba de la tendencia media-varianza encontrados por
##  modelGeneVarWithSpikes() por batches
hvg.pbmc.cv2.3 <- getTopHVGs(dec.cv2.pbmc, # Also works with modelGeneCV2()
                             var.field = "ratio", var.threshold = 1 # note `var.field` and value of `var.threshold`
)
str(hvg.pbmc.cv2.2)


### EJERCICIO: Dibujando los HVGs ###

## Repetir la gráfica que muestra la tendencia de la relación media-varianza
## (ejeX: media de la expresión, ejeY: varianza de la expresión) incluyendo la
## línea de tendencia (obtenida con modelGeneVar, modelGeneVarWithSpikes, modelGeneCV2).
## En esta gráfica, deberás colorear los puntos que corresponden a los HGVs obtenidos
## con algunos de los enfoques revisados


## Graficando los top 1000 genes con la mayor varianza encontrados por modelGeneVar()
## en verde
plot(fit.pbmc$mean, fit.pbmc$var,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression"
)
points(fit.pbmc$mean[hvg.pbmc.var], fit.pbmc$var[hvg.pbmc.var], col = "green")
curve(fit.pbmc$trend(x), col = "red", add = TRUE, lwd = 2)

## Graficando los genes con un p-value ajustado menor a 0.05 encontrados por modelGeneVar()
## en rosa
plot(fit.pbmc$mean, fit.pbmc$var,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression"
)
points(fit.pbmc$mean[hvg.pbmc.var.2], fit.pbmc$var[hvg.pbmc.var.2], col = "pink")
curve(fit.pbmc$trend(x), col = "blue", add = TRUE, lwd = 2)

## Graficando los genes por arriba de la tendencia media-varianza encontrados por modelGeneVar()
## en marrón
plot(fit.pbmc$mean, fit.pbmc$var,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression"
)
points(fit.pbmc$mean[hvg.pbmc.var.3], fit.pbmc$var[hvg.pbmc.var.3], col = "brown")
curve(fit.pbmc$trend(x), col = "orange", add = TRUE, lwd = 2)
