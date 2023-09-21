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

