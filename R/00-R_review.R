########## Breve repaso de R ##########

## Crear un data.frame
df <- data.frame(x = c(TRUE, FALSE, NA, NA), y = c(12, 34, 56, 78))
row.names(df) <- letters[1:4]

df

## Acceder a los nombres de las columnas
colnames(df)

## Acceder a los nombres de las filas
rownames(df)

## Obtener información booleana
df$y < 20
# Acceder al data.frame con la información booleana
df[df$y < 40, ]

## Acceder a información contenida dentro de los datos - %in% (dentro de)
bool_info <- rownames(df) %in% c("a", "c", "z")
df[bool_info, ]

## Acceder a información contenida dentro de los datos - & (y)
bool_info <- df$y < 50 & df$y > 20
df[bool_info, ]

## Acceder a información contenida dentro de los datos - | (o)
bool_info <- df$y < 20 | df$y > 60
df[bool_info, ]
