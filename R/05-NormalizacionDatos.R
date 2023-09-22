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
## de los clusters, solo la magnitud.

## La normalización por Library Size factor suele ser suficiente en algunas ocasiones
##donde se busca identificar clusters y los marcadores de los clusters.







