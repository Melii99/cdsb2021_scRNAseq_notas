########## Normalización de los Datos ##########

# La normalización tiene como objetivo remover las diferencias sistemáticas que
# surgen a partir de los procedimientos experimentales, para que no interfieran
# cuando comparamos los perfiles de expresión entre células.

# Al normalizar los datos, las diferencias observadas entre poblaciones célulares
# o condiciones son debido a la biología y no por factores técnicos.


### Conceptos Básicos ###

## Ejemplos de sesgos técnicos (TIP: ¿Qué es RPKM? - R. "Reads Per Kilobase Million):
# Sesgos de Composición, de Amplificación de PCR, de Posición en el Fragmento,
# de Longitud de Fragmento, de Captura de RNA...

## ¿Qué es correción por lote (batch effect correction)? Da un ejemplo.


## ¿Cuáles son las diferencias entre correción por lote y normalización?
# El objetivo del batch effect correction es mitigar las diferencias en condiciones
# experimentales (otra persona, hora o días diferentes, reactivos diferentes, etc)
# mientras que la normalización se refiere a "quitar" propiedades biofísicas
# (longitud del gen, profundidad de secuenciación, etc)


### Dataset zeisel ###

library("scRNAseq")
library("scater")

## Cargar los datos en R
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)
sce.zeisel

### QC ###

## Identificar los genes mitocondriales
is.mito <- which(rowData(sce.zeisel)$featureType == "mito")
stats <- perCellQCMetrics(sce.zeisel, subsets = list(Mt = is.mito))
qc <-
  quickPerCellQC(stats,
                 percent_subsets = c("altexps_ERCC_percent","subsets_Mt_percent")
                 )


sce.zeisel <- sce.zeisel[, !qc$discard]

## Cuántos genes son mitocondriales? (34)
length(is.mito)

## ¿Cuántos genes tienen: bajas cuentas de librería, bajos features,
## alto porcentaje de expresión de ERCC, alto porcentaje de genes MT?
## ¿Cuántas células descartamos? (TIP: perCellQCMetrics y quickPerCellQC)
colSums(as.data.frame(qc))


### scaling normalization ###

## La normalización por escalamiento es la estrategia más simple y usada.
## Representa el estimado del sesgo relativo en cada célula.
## Se realiza dividiendo todas las cuentas de cada célula por un factor de escalamiento
## específico para cada una. Factor de escalamiento - Library Size factor.
## (cuentas normalizadas = cuentas/ LibrarySizefactor)
## Suposición: Cualquier sesgo específico en cada célula (e.j. eficiencia en la
## captura o en la amplificación) afecta a todos los genes de igual manera a través de
## escalar por el promedio esperado de cuentas para dicha célula.


## Library Size: La suma total de las cuentas a tráves de todos los genes en una célula.
## sum(genes)

## Estimar tamaños de librerías
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)

## Examinar la distribución de los tamaños de librerías que acabamos de estimar
summary(lib.sf.zeisel)

## Histograma de la distribución de los tamaños de librerías
hist(log10(lib.sf.zeisel), xlab = "Log10[Library Size factor]", col = "grey80")


### Ejercicio: library Size ###


## Calcula library Size ls.zeisel
ls.zeisel <- colSums(counts(sce.zeisel))
summary(ls.zeisel)

## ¿Son idénticos ls.zeisel y lib.sf.zeisel?
## R. No, son "proporcionales" pero no iguales
identical(ls.zeisel, lib.sf.zeisel)

## ¿Son proporcionales?
plot(
  ls.zeisel,
  lib.sf.zeisel,
  log = "xy",
  main = "Proporcionalidad",
  xlab = "Library Size",
  ylab = "Library Size Factor"
)

## Calcula lib.sf.zeisel de forma manual
lib_size_factors <- ls.zeisel/mean(ls.zeisel)
summary(lib_size_factors)
identical(lib.sf.zeisel, lib_size_factors)

## Normalizar por Library Size factor asume que no hay desigualdad en la cantidad
## de genes differencialmente expresados (DE) entre dos células.
## Es decir, que para cada grupo de genes sobre-expresados, debe existir un grupo
## de genes sub-expresados en la misma magnitud, cuando esto no pasa se le conoce
## como sesgo de composición

## Para análisis exploratorios, la precisión de la normalización no es un punto
## mayor a considerar. El sesgo por composición normalmente no afecta la separación
## de los clusters, solo la magnitud (suele ser suficiente para la exploración).

## La normalización por Library Size factor suele ser suficiente en algunas ocasiones
##donde se busca identificar clusters y los marcadores de los clusters.



### Normalización por decircunvolución (deconvolution) ###

## Un sesgo técnico que es importante considerar es el sesgo
## de composición del transcriptoma (RNA)
## EJ. un gen X (o grupo de genes) se expresa en mayor cantidad en la célula A
## comparado a la célula B. Esto significa que más recursos fueron tomados por el
## gen X, disminuyendo la covertura de los demás

## En bulk RNA-seq Se assume que la mayoría de genes no estarán DE entre las muestras
##  y cualquier diferencia entre los genes non-DE representa un sesgo el cual se
## remueve (calculando un factor de normalización)

## En la Normalización por decircunvolución las células se

## Normalización por decircunvolución (deconvolution)
library("scran")

## Pre-clustering
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)

## Calcular factores de tamaño para la decircunvolución (deconvolution)
deconv.sf.zeisel <- calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)

## Examinar la distribución de los factores de tamaño
summary(deconv.sf.zeisel)

## Graficar los Deconvolution size factors
hist(log10(deconv.sf.zeisel),
     xlab = "Log10[Deconvolution size factor]",
     col = "grey80"
)

## Comparar los factores de normalización (size factors) obtenidos por
## escalamiento (lib.sf.zeisel) y deconvolution (deconv.sf.zeisel)
plot(lib.sf.zeisel,
     deconv.sf.zeisel,
     xlab = "Library size factor",
     ylab = "Deconvolution size factor",
     log = "xy",
     pch = 16
)
abline(a = 0, b = 1, col = "red")


### Ejercicios: deconvolution ###

## ¿Cúantos clusters rápidos obtuvimos? R. 12 clusters
levels(clust.zeisel)

## ¿Cúantas células por cluster obtuvimos?
## R. 113, 123, 140, 224, 231, 243, 252, 259, 281, 300, 324, 325
cells_cluster <- sort(table(clust.zeisel))
cells_cluster
barplot(cells_cluster)

## ¿Cúantos clusters rápidos obtendríamos si cambiamos el tamaño mínimo a 200?
## Usa 100 como la semilla (seed). R. 10 clusters
set.seed(100)
sort(table(quickCluster(sce.zeisel, min.size = 200)))

## ¿Cúantas líneas ves en la gráfica?
plot(lib.sf.zeisel,
     deconv.sf.zeisel,
     xlab = "Library size factor",
     ylab = "Deconvolution size factor",
     log = "xy",
     pch = 16,
     col = as.integer(factor(sce.zeisel$level1class))
)
abline(a = 0, b = 1, col = "red")
abline(a = -.2, b = 0.95, col = "red")
abline(a = 0.08, b = 1, col = "red")

## La normalización por decircunvolución (deconvolution) mejora los resultados
## para análisis posteriores de una manera más precisa que los métodos para bulk RNA-seq.

## scran algunas veces alcula factores negativos o ceros lo cual altera la matrix
## de expresión normalizada. ¡Checa los factores que calculas!
summary(deconv.sf.zeisel)

## Si obtienes factores negativos intenta variar el número de clusters, checa si
## incrementar el número de células por cluster te dan factores positivos.


### Transformación logarítmica ###

## En una escala logaritmica las diferencias entre la expresión de genes
## pueden verse mejor (más claramente comparables, ver diferencias sutiles)
## Ej. ¿Qué gen es más interesante?

## Gen X: el promedio de expresión en el tipo celular A: 50 y B: 10
50 - 10
log(50) - log(10)

## Gen Y: el promedio de expresión en el tipo celular A: 1100 y B: 1000
1100 - 1000
log(1100) - log(1000)


## Una vez calculados los factores de normalización con computeSumFactors(),
## podemos calular las cuentas en escala logaritmica usando logNormCounts().

## Los valores resultantes son valores de expresión normalizados transformados
## en escala logarítmica.

## Normalización de los datos
# Normalization
# set.seed(100)
# clust.zeisel <- quickCluster(sce.zeisel)
# sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clust.zeisel, min.mean=0.1)

## Log transformation (las cuentas normalizadas se guardan en assays como logcounts)
sce.zeisel <- scater::logNormCounts(sce.zeisel)
assayNames(sce.zeisel)


### Ejercicio: Transformación logarítmica ###

## ¿Qué es una pseudo-cuenta? R. Un número que se agrega para obtener el logaritmo
## ¿Porqué se usa? R. para evitar errores con log(0) (log(0) = -Inf)
## ¿Qué valor de pseudo-cuenta usa logNormCounts()? R. pseudo.count = 1


### Funciones interesantes para después de normalizar ###

## PCA de las cuentas de los datos de expresión génica del sce sce.zeisel
sce.zeisel <- runPCA(sce.zeisel)
## Graficar PC1 PC2
plotPCA(sce.zeisel, colour_by = "level1class")
## Graficar el  Relative Log Expression (RLE) (distribución de las diferencias en la expresión)
plotRLE(sce.zeisel, exprs_values = "logcounts", colour_by = "level1class")






