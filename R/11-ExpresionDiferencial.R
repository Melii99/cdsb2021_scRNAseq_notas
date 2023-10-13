########## Análisis de expresión diferencial ##########

## Paquetes usados ##
library("MouseGastrulationData") ## para descargar datos de ejemplo
library("scater") ## para gráficas y control de calidad
library("scran") ## para selección de genes, clustering, etc
library("batchelor") ## para métodos de correción de batch (lote)
library("patchwork") ## para agrupar gráficas
library("Polychrome") ## para muchos colores
library("bluster") ## métodos de clustering
library("edgeR") ## para expresión diferencial



## scRNA-seq nos puede ayudar a estudiar cambios en composición (cambios en
## proporciones de células) o cambios en niveles de expresión de genes entre
## varias condiciones biológicas


## El primero se llama cambios de abundancia,
# Ejemplo: después de un tratamiento con una droga
# Ejemplo: después de modificaciones genéticas

## Nos permite obtener mayor resolución biológica que experimentos convencionales
## de RNA-seq, sobre todo si podemos asociar cambios en poblaciones celulares a
## manipulaciones experimentales


### Dos categorías de análisis ###


## Análisis de expresión diferencial ##

## Buscamos cambios en niveles de expresión entre condiciones para células del
## mismo tipo que están presentes en todas las condiciones !!!


## Análisis de abundancia diferencial ##

## Buscamos cambios en la composición de los tipos celulares entre condiciones !!!

## Podría ser entre estados celulares en vez de tipos celulares !!!

## Son dos lados de la misma moneda



### Datos de ejemplo ###

## Embriones de ratón quiméricos. Pijuan-Sala, B. et al. A single-cell molecular
## map of mouse gastrulation and early organogenesis. Nature 566, 490–495 (2019).

## 3 réplicas: batches


## Descarga de los datos desde bioconductor y guardar el SCE
library("MouseGastrulationData")
sce.chimera <- WTChimeraData(samples = 5:10)

## Explorar los datos: resumen estadístico para cada columna
sapply(colData(sce.chimera)[, -(1:2)], function(x) {
  x <- if (is.character(x) || is.integer(x)) factor(x) else x
  summary(x)
})

## Encontramos que:

## sample: 6 ratones diferentes
## tomato: inyectados o no con td-Tomato
## pool: lote de secuenciación, cada lote con 1 con y otro sin inyección
## celltype.mappped: 35 tipos de células anotados


## Número de células en nuestras variables principales ##
## (tabla con la distribución de células en función de las combinaciones únicas
## de los valores de sample, pool y tomato)
with(colData(sce.chimera), table(sample, pool, tomato))


## Número de tipos celulares ##
## Cálculo del número de tipos celulares únicos presentes en sce.chimera
length(unique(sce.chimera$celltype.mapped))



### Procesamiento de los datos ###

## Usaremos batchelor porque tenemos muestras de 3 lotes de muestras y
## queremos eliminar diferencias entre los lotes !!!


## Anotacion ##

library("scater")
rownames(sce.chimera) <- uniquifyFeatureNames(
  rowData(sce.chimera)$ENSEMBL, rowData(sce.chimera)$SYMBOL
)


## Control de calidad QC ##

## Filtrar y eliminar las células ya marcadas como "stripped" o "Doublet"
drop <- sce.chimera$celltype.mapped %in% c("stripped", "Doublet")
sce.chimera <- sce.chimera[, !drop]


## Normalizacion ##
sce.chimera <- logNormCounts(sce.chimera)


## Modelado de la varianza ##

## Modelar la varianza yseleccionar los genes con variabilidad biológica distinta de 0 (HVGs)
library("scran")
dec.chimera <- modelGeneVar(sce.chimera, block = sce.chimera$sample)
chosen.hvgs <- dec.chimera$bio > 0


## ¡¡¡ Merging !!! ##

## Uso de batchelor ##

## corregir batch effects y combinar los datos de diferentes experimentos
library("batchelor")
set.seed(01001001)
merged <- correctExperiments(sce.chimera, #SCE
                             batch = sce.chimera$sample, # batches
                             subset.row = chosen.hvgs, # HVGs
                             PARAM = FastMnnParam( # Orden de los datos
                               merge.order = list( # Orden en que se combinan
                                 list(1, 3, 5), # WT (3 replicates)
                                 list(2, 4, 6) # td-Tomato (3 replicates)
                               )
                             )
)


## Clustering ##

## Creación de un grafo SNN
g <- buildSNNGraph(merged, use.dimred = "corrected")

## Encontrar los clusters aplicando el algoritmo de Louvain
clusters <- igraph::cluster_louvain(g)

## Asignación de etiquetas de colores a las células (merged) con su respectivo clúster
colLabels(merged) <- factor(clusters$membership)


## Reducción de dimensionalidad ##

## Reducción dimensionalidad tSNE
merged <- runTSNE(merged, dimred = "corrected", external_neighbors = TRUE)

## Reducción dimensionalidad UMAP
merged <- runUMAP(merged, dimred = "corrected", external_neighbors = TRUE)



### Exploración de los datos de ejemplo ###


## Exploremos si tenemos clusters con una diferencia grande en el número de celulas
## entre las muestras sin y con inyecciones de td-Tomato

## Exploremos el número de células en cada cluster a lo largo de los
## 3 lotes de secuenciación (batch)


## Clusters vs DE por td-Tomato
table(colLabels(merged), merged$tomato)

## Clusters vs lotes de muestras (batch)
table(colLabels(merged), merged$pool)

## Visualización ##

## Visualizar nuestros clusters que (26 en dimensiones reducidas) de t-SNE
## Queremos que todos los clusters tengan muestras de cada lote de secuenciación (batch).
## Vemos que no parece que haya mucha señal en base a td-Tomato
library("patchwork")
plotTSNE(merged, colour_by = "tomato", text_by = "label") +
  plotTSNE(merged, colour_by = data.frame(pool = factor(merged$pool)))

## Podemos usar facet_wrap() para reducir el over-plotting y ver mejor la información:
##  t-SNE donde las células se agrupan por las categorías de tomato y se visualizan en
## paneles separados para cada categoría
plotTSNE(merged,
         colour_by = "tomato",
         other_fields = c("tomato")
) + facet_wrap(~tomato)

##  t-SNE donde las células se agrupan por las categorías de pool y se visualizan en
## paneles separados para cada categoría
plotTSNE(merged,
         colour_by = data.frame(pool = factor(merged$pool)),
         other_fields = c("pool")
) + facet_wrap(~pool)



### Nuestros clusters vs los originales ###


## Las siguientes gráficas nos ayudan a comparar nuestros clusters vs los que
## encontraron en el estudio original

## Para cols_label
## Definir colores, de no hacerlo scater nos los pone en una escala continua
cols_label <- Polychrome::palette36.colors(length(unique(merged$label)))
## Asignar las labels a los colores correspondientes en la paleta
names(cols_label) <- unique(merged$label)

# Para cols_celltype.mapped
## Definir colores, de no hacerlo scater nos los pone en una escala continua
cols_celltype.mapped <- Polychrome::palette36.colors(length(unique(merged$celltype.mapped)))
## Asignar las labels a los colores correspondientes en la paleta
names(cols_celltype.mapped) <- unique(merged$celltype.mapped)


## Tsne de Nuestros clusters encontrados vs anotación de células por los autores originales
plotTSNE(merged, colour_by = "label", text_by = "label") +
  theme(legend.position = "none") +
  scale_colour_manual(values = cols_label) +
  plotTSNE(merged, colour_by = "celltype.mapped") +
  theme(legend.position = "none") +
  scale_colour_manual(values = cols_celltype.mapped)

## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.
## Scale for 'colour' is already present. Adding another scale for 'colour',
## which will replace the existing scale.


## ¿Parecen similares?

## Es difícil el proceso de comparar clusters (a simple vista)

## Podemos usar bluster para evaluar númericamente que tanto se parecen los clusters.
## Entre más cerca de 1, mejor en pairwiseRand()

## También podemos hacer un heatmap


## Evaluar numéricamente el parecido entre los clusters
library("bluster")
pairwiseRand(colLabels(merged), merged$celltype.mapped, "index")

## Heatmap de las coincidencias entre los clusters encontrados y los originales
by.label <- table(colLabels(merged), merged$celltype.mapped)
pheatmap::pheatmap(log2(by.label + 1), color = viridis::viridis(101)) # Log2 y +1 para manejar valores cero



### Análisis de expresión diferencial ###


## En RNA-seq estamos acostumbrados a evaluar si hay diferencias en los niveles  de expresión
## de genes entre condiciones, así que es natural que lo hagamos con scRNA-seq también

## Pero los datos de scRNA-seq tienen muchos ceros



## Pseudo-bulking ##

## El proceso de pseudo-bulking es un truco que nos permite usar métodos de bulk
## RNA-seq para analizar nuestros datos de scRNA-seq

## Se suma la expresión de diferentes células por mismo gen y misma condición
## (eliminar problema de 0s  comprimir la matriz)

## Cómo tenemos muchas células de cada condición, para cada gene podemos sumar los
## niveles de expresión entre todas las células de esa condición


## Ejemplo de mi trabajo:
  # 12 muestras
  # 7 regiones
  # 47,681 spots (digamos que células)

# Podemos comprimir la información a una matriz de 12 * 7 = 84 columnas

# Nos quedamos con pocas réplicas para nuestro análisis, pero justamente los métodos
# de bulk RNA-seq están diseñados para esos escenarios (claro, entre más datos mejor!!!)


## Podemos hacerlo manualmente o de forma más sencilla con la función aggregateAcrossCells() ##

## La función aggregateAcrossCells() agrupa los datos de merged
## según las combinaciones únicas de "celltype.mapped" y "sample" y realiza una
## agregación (como sumar o promediar) dentro de cada grupo. Los resultados de
## esta operación de agregación se almacenan en el objeto summed. Dependiendo del contexto,
## esta agregación puede ser útil para resumir los datos y facilitar análisis posteriores

# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
summed <- aggregateAcrossCells(merged,
                               id = colData(merged)[, c("celltype.mapped", "sample")]
)

summed

## Dimensiones iniciales
dim(merged)

## Dimensiones reducidas
dim(summed)

## Numero teorico de columnas final, pero muchas columnas npo tienen  datos, por eso 186
with(colData(merged), length(unique(celltype.mapped)) * length(unique(sample)))

## En teoría podríamos tener más columnas, pero no las tenemos todas porque no
## tenemos datos para todas las combinaciones

## Esto puede afectar nuestro análisis, y pues afecta cuantas variables podremos
## usar para ajustar

## Por ejemplo, si agregamos sexo con 2 opciones, duplicaríamos el número teórico
## de columnas pero tal vez no tengamos suficientes datos

## Si lo llevas al extremo, terminas con los mismos datos de scRNA-seq que con los que empezamos



## Convertir a un objeto nuevo ##

## Hagamos nuestro análisis de expresión diferencial

## Empezaremos con solo un tipo celular: Mesenchyme

label <- "Mesenchyme"
current <- summed[, label == summed$celltype.mapped]
dim(current)


## Vemos que nos quedamos con solo 14,699 genes a lo largo de 6 muestras
## Esto sería un experimento pequeño de bulk RNA-seq



### Edge R ###

## Usaremos edgeR de Robinson, McCarthy e Smyth, Bioinformatics, 2010 que es uno
## de los paquetes más usados para análisis de expresión diferencial en bulk RNA-seq

## Crear una DGEList para usar en EdgeR
library("edgeR")
y <- DGEList(counts(current), samples = colData(current))
y


## Pre-procesamiento ##

## Antes de poder continuar, vamos a eliminar muestras que construimos con el
## proceso de pseudo-bulking que no tengan al menos 10 células
discarded <- current$ncells < 10
y <- y[, !discarded]
summary(discarded)

## Eliminar los genes que tengan bajos niveles de expresión
keep <- filterByExpr(y, group = current$tomato)
y <- y[keep, ]
summary(keep)


## Normalizar los datos ##

## Pero si ya habíamos normalizado los datos de scRNA-seq, ¿qué pasó?
## Empezamos de nuevo con las cuentas originales en y <- DGEList(counts(current),
## samples=colData(current)) y no las normalizadas.

## calcNormFactors() asume que la mayoría de los genes no están diferencialmente expresados
y <- calcNormFactors(y)
y$samples

## Podemos visualizar los cambios de expresión para todos los genes, una muestra a la vez
par(mfrow = c(2, 3))
for (i in seq_len(ncol(y))) {
  plotMD(y, column = i)
}

## Podemos visualizar los cambios de expresión para todos los genes, una muestra a la vez
par(mfrow = c(2, 3))
for (i in seq_len(ncol(y))) {
  plotMD(y, column = i)
}

## Podemos repetir el plotMDS() pero con colores por lote (batch) de pool.
plotMDS(cpm(y, log = TRUE),
        col = c("3" = "darkorchid1", "4" = "darkblue", "5" = "tomato4")[factor(y$samples$pool)]
)

## Acá vemos que si hay diferencias entre lotes, en particular entre el lote de las
## muestras 1 y 2 y el resto, ya que el eje X explica el 38% de la varianza.



### Modelo estadístico ###


## Si todo nos parece bien, podemos seguir con definir nuestro modelo estadístico

## Vamos a ajustar por lote (batch) y encontrar diferencias por la inyección de td-Tomato

## Como empezamos con las cuentas desde cero, tenemos que tomar en cuenta la
## variación por lote de secuenciación

design <- model.matrix(~ factor(pool) + factor(tomato), y$samples)
design

## Si queremos explorar nuestro modelo estadístico de forma interactiva,
## podemos usar ExploreModelMatrix
if (interactive()) {
  ExploreModelMatrix::ExploreModelMatrix(y$samples[, c("pool", "tomato")], ~ factor(pool) + factor(tomato))
}

## Tal y como en bulk RNA-seq, podemos usar la información de los genes para
## mejorar nuestros estimados de la varianza para cada gene, de tal forma que
## mejoramos los resultados estadísticos aunque tengamos pocas muestras
y <- estimateDisp(y, design)
summary(y$trended.dispersion)

## Visualizando la varianza (dispersión estimada)
plotBCV(y)

## Ajustar el modelo lineal generalizado y realiza una estimación robusta de los parámetros
fit <- glmQLFit(y, design, robust = TRUE)

## Resumen de la estimación de la varianza del modelo
summary(fit$var.prior)

## Resumen estadístico de sobre los grados de libertad previos del modelo
summary(fit$df.prior)

## Evaluar la calidad del ajuste del modelo en términos de dispersión
plotQLDisp(fit)

## Modelo estadistico ##

## Ahora si podemos correr nuestro modelo estádistico
## Se ejecutan pruebas de contraste de tipo LRT (Likelihood Ratio Test) para
## identificar genes diferencialmente expresados en un modelo ajustado utilizando
## el paquete edgeR en R.
res <- glmQLFTest(fit, coef = ncol(design))
de_n <- summary(decideTests(res))
de_n

## Resumen de los genes encontrados como diferencialmente expresados
topTags(res)

## Encontramos 16 genes diferencialmente expresados por la inyección de td-Tomato.



### Análisis de expresión diferencial - De forma sencilla ###


## La función pseudoBulkDGE() corre todos esos pasos por nosotros!!!

## Remover todas las muestras de pseudo-bulk con un numero 'insuficiente' de células
summed.filt <- summed[, summed$ncells >= 10]

## Realizar un análisis de expresión génica diferencial en datos pseudo-bulk
library("scran")
de.results <- pseudoBulkDGE(summed.filt, # objeto
                            label = summed.filt$celltype.mapped, # Grupos
                            design = ~ factor(pool) + tomato, # Relación entre variables predictoras (modelo lineal)
                            coef = "tomatoTRUE", # Contraste a evaluar
                            condition = summed.filt$tomato # Condiciones exp.
)

## Tipo de dato de la salida (Nos regresa una lista con los resultados para cada uno de nuestros tipos celulares)
class(de.results)

length(de.results)

## Podemos extraer los resultados para nuestro tipo celular de interés !!!
## por ejemplo Allantois.

## Guardar resultados de expresión génica diferencial específicamente para el tejido "Allantois"
cur.results <- de.results[["Allantois"]]
## Ordenar estos resultados por p-value (ascendente)
cur.results[order(cur.results$PValue), ]

## Guardar el coeficiente biológico de variación específico para "Allantois"
y.allantois <- metadata(cur.results)$y
## Visualizarlo
plotBCV(y.allantois)


## También nos dice que tipos celulares fallaron porque no teníamos suficiente
## información para hacer el análisis
metadata(de.results)$failed

## Podemos hacer la misma gráfica que hicimos de forma manual para Mesenchyme.

## Guardar resultados de expresión génica diferencial específicamente para el tejido "Mesenchyme"
cur.results.Mesenchyme <- de.results[["Mesenchyme"]]
## Guardar el coeficiente biológico de variación específico para "Mesenchyme"
y.Mesenchyme <- metadata(cur.results.Mesenchyme)$y
## Visualizar
plotBCV(y.Mesenchyme)
