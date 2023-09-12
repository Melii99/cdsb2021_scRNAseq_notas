########## Ejercicio 2 ##########

## Importamos librerías
library("SingleCellExperiment")
library("scRNAseq")

## Mini muestreo del set de datos usado en:
#https://bioconductor.org/books/release/OSCA/zeisel-mouse-brain-strt-seq.html#introduction-5

## Descarga de archivos
archivo_cuentas <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/min_sce.csv"
archivo_rowData <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/rowD.csv"
archivo_colData <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/colD.csv"

## Lectura de los datos de los archivos descargados
counts <- read.csv(archivo_cuentas, row.names = 1, header = TRUE, check.names = F)
col.data <- DataFrame(read.csv(archivo_colData, row.names = 1, header = TRUE, check.names = F))
row.data <- read.csv(archivo_rowData, row.names = 1, header = TRUE, check.names = F)

## Crear un objeto SingleCellExperiment
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = col.data,
  rowData = row.data
)

## Agregar logcounts
sce <- scater::logNormCounts(sce)

sce # 100 columnas-células, 30 filas-genes

## ¿Cuántos genes tenemos? (30)
ncol(assay(sce))

## ¿Qué información tenemos en el rowData()?

# Nombres de las columnas de rowData()
colnames(rowData(sce))
# Vistazo de los datos de rowData()
head(rowData(sce))

## Genes de interés
int_gen <- c("Angpt1", "Chic2", "Mir503", "Magee2", "Nenf", "Eps15l1", "Hsf2bp", "Gnptg", "Vegfb", "Atmin", "Gad1", "Gad2", "Slc32a1", "Dner", "Slc2a13", "Slc6a1", "Nrxn3")

## Extraer los datos de los genes que nos interesan (objeto int_gen)

# rowData(sce): Accede a la información de los genes (nombres de genes, símbolos, anotaciones, etc.)
#rownames(): Accede a los nombres de las filas (en rowData)
# %in%: Identificadores de genes en rowData que están presentes en el vector int_gen (genes de interés)
# [filas, columnas]

sce[rownames(rowData(sce)) %in% int_gen, ]


## Crear un objeto llamado min_sce con los datos de solo esos genes
min_sce <- sce[rownames(rowData(sce)) %in% int_gen, ]

## ¿Cuáles son parte del tejido interneurons o del tejido pyramidal CA1 ? (del objeto min_sce)

# Observamos la informacion de colData (cell metadata)
colData(sce) # descubrimos que "interneurons" y "pyramidal CA1" se encuentran en la columna level1class

min_sce[, min_sce$level1class == "interneurons" | min_sce$level1class == "pyramidal CA1"]

## Con este subconjunto, crea el objeto tej_min_sce
tej_min_sce <- min_sce[, min_sce$level1class == "interneurons" | min_sce$level1class == "pyramidal CA1"]

## Una vez que tengan el objeto ´SingleCellExperiment´ llamado ´tej_min_sce´, corran el siguiente código.

library("scater")
plotHeatmap(object = tej_min_sce, features = rownames(tej_min_sce), order_columns_by = "level1class")


## Ejercicio extra

## Realizar los mismos pasos, pero ahora con el Ensembl gene IDs de la lista de genes

# Lista de Ensembl gene IDs de interés
int_gen <- c("ENSMUSG00000071076", "ENSMUSG00000002076", "ENSMUSG00000024962", "ENSMUSG00000031224", "ENSMSG00000036560", "ENSMUSG00000037499", "ENSMUSG00000006276", "ENSMUSG00000035521", "ENSMUSG00000047388", "ENSMUSG0000051079", "ENSMUSG00000076122", "ENSMUSG00000029229", "ENSMUSG00000022309", "ENSMUSG00000036766", "ENSMUSG00000070880", "ENSMUSG00000026787", "ENSMUSG00000066392", "ENSMUSG00000036298", "ENSMUSG00000037771", "ENSMUSG00000030310")

# Extracción de datos de genes de interés utilizando los Ensembl gene IDs
sce[rowData(sce)$Ensembl %in% int_gen, ]

# Crear un objeto llamado min_sce con los datos de solo esos genes (usando rowData(sce)$Ensembl para ensembl IDs)
min_sce <- sce[rowData(sce)$Ensembl %in% int_gen, ]

# Identificación de genes que pertenecen al tejido interneurons o al tejido pyramidal CA1
min_sce[, min_sce$level1class == "interneurons" | min_sce$level1class == "pyramidal CA1"]

# Crear un objeto llamado tej_min_sce con las células seleccionadas
tej_min_sce <- min_sce[, min_sce$level1class == "interneurons" | min_sce$level1class == "pyramidal CA1"]

## Una vez que tengan el objeto ´SingleCellExperiment´ llamado ´tej_min_sce´, corran el siguiente código.

library("scater")
plotHeatmap(object = tej_min_sce, features = rownames(tej_min_sce), order_columns_by = "level1class")


### Información de la sesión de R ###
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

