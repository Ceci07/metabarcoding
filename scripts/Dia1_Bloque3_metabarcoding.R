# ============================================================
# TALLER DE ANÁLISIS DE DATOS DE METABARCODING
# Día 1 - Bloque 3: conocer los datos del taller
# Caso de estudio: Zwirzitz et al. (2020)
# ============================================================
#
# En este bloque vamos a:
#   1. Descargar la metadata y los archivos FASTQ desde Google Drive.
#   2. Explorar la información de las muestras.
#   3. Reconocer los archivos forward (_1) y reverse (_2).
#
# Todavía NO vamos a usar DADA2, mia ni phyloseq.


# ------------------------------------------------------------
# 1. Revisar el entorno de trabajo
# ------------------------------------------------------------

R.version.string
getwd()

# Todo lo que descarguemos quedará dentro de la carpeta "data".
dir.create("data", showWarnings = FALSE)

cat("\nPREGUNTA: ¿En qué directorio se creó la carpeta data?\n")


# ------------------------------------------------------------
# 2. Conectarse a la carpeta pública de Google Drive
# ------------------------------------------------------------

if (!requireNamespace("googledrive", quietly = TRUE)) {
  install.packages("googledrive")
}

library(googledrive)

 # No hace falta iniciar sesión porque la carpeta es pública.
drive_deauth()


carpeta_taller <- as_id(
  "https://drive.google.com/drive/folders/1SOhCrJXgrtucWavxsxizLY4FErL3Almk"
)

# Vemos qué contiene la carpeta principal.
archivos_principales <- drive_ls(carpeta_taller)
archivos_principales


# ------------------------------------------------------------
# 3. Descargar y leer la metadata
# ------------------------------------------------------------

metadata_local <- file.path("data", "sampleData_alimentos.csv")

# Descargamos la metadata solamente si todavía no está en la computadora.
if (!file.exists(metadata_local)) {
  metadata_drive <- archivos_principales[
    archivos_principales$name == "sampleData_alimentos.csv",
  ]

  if (nrow(metadata_drive) == 0) {
    stop("No se encontró sampleData_alimentos.csv en Google Drive.")
  }

  drive_download(
    metadata_drive,
    path = metadata_local,
    overwrite = FALSE
  )
}

# Leemos la tabla. 

metadata <- read.csv(
  metadata_local,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# ------------------------------------------------------------
# 4. Explorar la metadata
# ------------------------------------------------------------

head(metadata)

# Cada fila corresponde a una muestra. Las columnas son:
#   Run: identificador del archivo en el repositorio de secuencias.
#   Sample: nombre de la muestra usado en el estudio.
#   Sampling position: lugar o etapa del proceso donde se tomó.
#   Farmer: productor de origen; vacío para muestras ambientales.



cat("\nEJERCICIO 1: ¿Cuántas muestras y cuántas variables hay?\n")
dim(metadata)
names(metadata)
str(metadata)


cat("EJERCICIO 2: ¿Cuántas posiciones de muestreo diferentes aparecen?\n")
table(metadata$`Sampling.position`)
length(table(metadata$Sampling.position))

cat("EJERCICIO 3: ¿Por qué algunas muestras no tienen productor (Farmer)?\n")
table(metadata$Farmer, useNA = "ifany")

table(
  metadata$Sampling.position[is.na(metadata$Farmer)],
  useNA = "ifany"
)

       
# ------------------------------------------------------------
# 5. Entrar a la subcarpeta "alimentos" y listar los FASTQ
# ------------------------------------------------------------

# Los FASTQ no están en la carpeta principal: están dentro de "alimentos".
carpeta_alimentos <- archivos_principales[
  archivos_principales$name == "alimentos",
]

if (nrow(carpeta_alimentos) == 0) {
  stop("No se encontró la subcarpeta 'alimentos' en Google Drive.")
}

archivos_alimentos <- drive_ls(as_id(carpeta_alimentos$id[1]))

# Nos quedamos solamente con archivos .fastq, .fq y sus versiones comprimidas.
es_fastq <- grepl(
  "\\.(fastq|fq)(\\.gz)?$",
  archivos_alimentos$name,
  ignore.case = TRUE
)

fastq_drive <- archivos_alimentos[es_fastq, ]

cat("\nCantidad de FASTQ encontrados:", nrow(fastq_drive), "\n")
print(fastq_drive$name)


# ------------------------------------------------------------
# 6. Separar forward y reverse y verificar los pares
# ------------------------------------------------------------

# En estos datos, _1 identifica las lecturas forward y _2 las reverse.
forward <- fastq_drive[grepl("_1\\.(fastq|fq)(\\.gz)?$",
                             fastq_drive$name, ignore.case = TRUE), ]

reverse <- fastq_drive[grepl("_2\\.(fastq|fq)(\\.gz)?$",
                             fastq_drive$name, ignore.case = TRUE), ]

cat("Archivos forward:", nrow(forward), "\n")
cat("Archivos reverse:", nrow(reverse), "\n")

# Quitamos _1.fastq.gz o _2.fastq.gz para obtener el ID de cada muestra.
id_forward <- sub("_1\\.(fastq|fq)(\\.gz)?$", "", forward$name,
                  ignore.case = TRUE)
id_reverse <- sub("_2\\.(fastq|fq)(\\.gz)?$", "", reverse$name,
                  ignore.case = TRUE)

sin_reverse <- setdiff(id_forward, id_reverse)
sin_forward <- setdiff(id_reverse, id_forward)

if (length(sin_reverse) == 0 && length(sin_forward) == 0 &&
    nrow(forward) == nrow(reverse) && nrow(forward) > 0) {
  cat("Todos los archivos forward tienen su correspondiente reverse.\n")
} else {
  cat("Hay pares incompletos. Revisar las siguientes muestras:\n")
  print(unique(c(sin_reverse, sin_forward)))
}

cat("\nPREGUNTA: ¿Por qué cada muestra tiene dos archivos FASTQ?\n")


# ------------------------------------------------------------
# 7. Descargar los FASTQ
# ------------------------------------------------------------

# Esta descarga deja preparados los datos para los bloques siguientes.
# Los archivos se guardan en data/fastq y no se vuelven a descargar si
# ya existen. Si solo queremos mirar la lista, cambiamos TRUE por FALSE.
descargar_fastq <- FALSE

if (descargar_fastq) {
  dir.create(file.path("data", "fastq"), showWarnings = FALSE)

  for (i in seq_len(nrow(fastq_drive))) {
    destino <- file.path("data", "fastq", fastq_drive$name[i])

    if (!file.exists(destino)) {
      cat("Descargando", fastq_drive$name[i], "\n")
      drive_download(
        fastq_drive[i, ],
        path = destino,
        overwrite = FALSE
      )
    } else {
      cat("Ya existe, no se descarga nuevamente:",
          fastq_drive$name[i], "\n")
    }
  }
}


# ------------------------------------------------------------
# 8. ¿Qué haremos después con estos datos?
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FASTQ\n")
cat("  -> DADA2: control de calidad y corrección de errores\n")
cat("  -> ASVs: secuencias observadas en las muestras\n")
cat("  -> Taxonomía: identificación de los microorganismos\n")
cat("  -> Diversidad: comparación de las comunidades microbianas\n")
cat("============================================================\n")


