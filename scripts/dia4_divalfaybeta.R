###############################################################################
# TALLER DE METABARCODING - DIA 4
# Diversidad alfa y diversidad beta
###############################################################################

# En este práctico compararemos dos sitios de muestreo:
#   - Sticking
#   - After singeing
#
# El objeto tse_dia4_base.rds fue preparado al finalizar el Día 3. Contiene:
# conteos, conteos rarefactados, abundancias relativas y presencia/ausencia.


# 1. PREPARAR LA SESION -------------------------------------------------------

# Seleccionar como directorio de trabajo la carpeta que contiene los datos.
# En RStudio también se puede usar: Session > Set Working Directory.

library(mia)
library(vegan)
library(ggplot2)

# Cargar el objeto tse de la carpeta de drive Día 4.
tse <- readRDS("tse_dia4_base.rds")

# Ver qué contiene el objeto.
tse

dim(assay(tse, "counts"))  
# Cantidad de ASVs (filas) y cantidad de muestras (columnas)

rownames(assay(tse, "counts"))[1:6]  
# Muestra los nombres de las primeras seis ASVs

colnames(assay(tse, "counts"))[1:6]  
# Muestra los identificadores de las primeras seis muestras

assayNames(tse)  # Muestra los nombres de las distintas tablas almacenadas en el objeto tse

assay(tse, "counts")[1:6, 1:6]  
# Conteos originales: número de lecturas de cada ASV en cada muestra

assay(tse, "relabundance")[1:6, 1:6]  
# Abundancias relativas: proporción que representa cada ASV dentro de cada muestra

assay(tse, "subsampled")[1:6, 1:6]  
# Conteos luego de la rarefacción: todas las muestras fueron igualadas a la misma profundidad

assay(tse, "rarefiedRelabundance")[1:6, 1:6]  
# Abundancias relativas calculadas a partir de los conteos rarefactados

assay(tse, "pa")[1:6, 1:6]  
# Presencia/ausencia: 1 indica que la ASV está presente y 0 que está ausente

# Un pequeño chequeo

colSums(assay(tse, "counts"))[1:6]  
# Cantidad total de lecturas originales en cada una de las primeras seis muestras

colSums(assay(tse, "relabundance"))[1:6]  
# Cada columna debería sumar 1 porque son proporciones

colSums(assay(tse, "subsampled"))[1:6]  
# Todas deberían sumar 6513: profundidad utilizada para rarefaccionar

colSums(assay(tse, "rarefiedRelabundance"))[1:6]  
# Cada muestra debería sumar 1 porque son abundancias relativas

table(assay(tse, "pa"))  
# Cuenta cuántos ceros (ausencias) y unos (presencias) hay en toda la tabla

# 2. SELECCIONAR LAS MUESTRAS -------------------------------------------------

# Para este ejercicio trabajaremos solamente con muestras de Sticking y de
# After singeing
tse_d4 <- tse[
  ,
  tse$Sampling.position %in% c("Sticking", "After singeing")
]
# Selecciona muestras de la superficie antes y después del chamuscado

# Eliminar categorías que ya no están presentes después de seleccionar muestras.
tse_d4$Sampling.position <- droplevels(tse_d4$Sampling.position)
tse_d4$Farmer <- droplevels(tse_d4$Farmer)

# Comprobar cuántas muestras quedaron en cada grupo.
table(tse_d4$Sampling.position)


# 3. DIVERSIDAD ALFA ----------------------------------------------------------

# La diversidad alfa describe la diversidad dentro de cada muestra.

# Riqueza observada: número de ASV detectadas.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "subsampled",
  index = "observed_richness",
  name = "observed"
)

# Shannon: combina riqueza y distribución de abundancias.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "subsampled",
  index = "shannon_diversity",
  name = "shannon"
)

# Uniformidad de Simpson: indica cuán uniformes son las abundancias.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "subsampled",
  index = "simpson_evenness",
  name = "simpson_evenness"
)

# Chao1: estima la riqueza considerando las ASV poco frecuentes.
tse_d4 <- addAlpha(
  tse_d4,
  assay.type = "counts",
  index = "chao1_richness",
  name = "chao1"
)

# Construir una tabla con los resultados.
alpha <- as.data.frame(colData(tse_d4))[
  ,
  c(
    "Sample", "Sampling.position", "observed", "shannon",
    "simpson_evenness", "chao1"
  )
]

alpha

# 4. COMPARAR LA DIVERSIDAD ALFA ---------------------------------------------

# Las muestras de Sticking y After singeing corresponden a los mismos
# cuatro animales. Por lo tanto, las observaciones están pareadas.

# Extraer el número que identifica al animal desde el nombre de la muestra.

alpha$animal <- gsub("\\D", "", alpha$Sample)


# Ordenar las etapas para mostrarlas en su orden temporal.

alpha$Sampling.position <- factor(
  alpha$Sampling.position,
  levels = c("Sticking", "After singeing")
)


# Comprobar que los animales estén en el mismo orden en ambas etapas.

animales_sticking <- alpha$animal[
  alpha$Sampling.position == "Sticking"
]

animales_singeing <- alpha$animal[
  alpha$Sampling.position == "After singeing"
]

pares_animales <- data.frame(
  Sticking = animales_sticking,
  After_singeing = animales_singeing
)

pares_animales


# Si los identificadores no coinciden, R detiene el análisis para evitar
# comparar entre sí muestras pertenecientes a animales diferentes.

stopifnot(
  animales_sticking == animales_singeing
)


# EJEMPLO: ¿qué haríamos si las muestras fueran independientes?
#
# En ese caso podríamos utilizar la fórmula:
#
# wilcox.test(
#   shannon ~ Sampling.position,
#   data = alpha,
#   exact = FALSE
# )
#
# Este comando no se utiliza aquí porque se siguieron los mismos
# animales antes y después del chamuscado.


# Aplicar la prueba de Wilcoxon pareada para cada índice.
# La prueba evalúa si existe un cambio sistemático dentro de los animales
# entre Sticking y After singeing.

wilcox_observed <- wilcox.test(
  alpha$observed[alpha$Sampling.position == "Sticking"],
  alpha$observed[alpha$Sampling.position == "After singeing"],
  paired = TRUE,
  exact = TRUE
)

wilcox_shannon <- wilcox.test(
  alpha$shannon[alpha$Sampling.position == "Sticking"],
  alpha$shannon[alpha$Sampling.position == "After singeing"],
  paired = TRUE,
  exact = TRUE
)

wilcox_simpson <- wilcox.test(
  alpha$simpson_evenness[alpha$Sampling.position == "Sticking"],
  alpha$simpson_evenness[alpha$Sampling.position == "After singeing"],
  paired = TRUE,
  exact = TRUE
)

wilcox_chao1 <- wilcox.test(
  alpha$chao1[alpha$Sampling.position == "Sticking"],
  alpha$chao1[alpha$Sampling.position == "After singeing"],
  paired = TRUE,
  exact = TRUE
)


# Mostrar los resultados individuales.

wilcox_observed
wilcox_shannon
wilcox_simpson
wilcox_chao1


# Reunir los valores p en una sola tabla.

resultados_wilcoxon <- data.frame(
  indice = c(
    "Riqueza observada",
    "Shannon",
    "Uniformidad de Simpson",
    "Chao1"
  ),
  p_value = c(
    wilcox_observed$p.value,
    wilcox_shannon$p.value,
    wilcox_simpson$p.value,
    wilcox_chao1$p.value
  )
)

resultados_wilcoxon


# Como analizamos cuatro índices, ajustamos los valores p mediante
# el procedimiento de Benjamini-Hochberg.

resultados_wilcoxon$p_ajustado_BH <- p.adjust(
  resultados_wilcoxon$p_value,
  method = "BH"
)

resultados_wilcoxon


# 4.1. ESTADÍSTICAS DESCRIPTIVAS ---------------------------------------------

# Media y desviación estándar: resumen basado en el promedio.

media_alpha <- aggregate(
  cbind(
    observed,
    shannon,
    simpson_evenness,
    chao1
  ) ~ Sampling.position,
  data = alpha,
  FUN = mean
)

sd_alpha <- aggregate(
  cbind(
    observed,
    shannon,
    simpson_evenness,
    chao1
  ) ~ Sampling.position,
  data = alpha,
  FUN = sd
)

media_alpha
sd_alpha


# Mediana y rango intercuartílico: resumen robusto y coherente
# con la comparación no paramétrica.

mediana_alpha <- aggregate(
  cbind(
    observed,
    shannon,
    simpson_evenness,
    chao1
  ) ~ Sampling.position,
  data = alpha,
  FUN = median
)

iqr_alpha <- aggregate(
  cbind(
    observed,
    shannon,
    simpson_evenness,
    chao1
  ) ~ Sampling.position,
  data = alpha,
  FUN = IQR
)

mediana_alpha
iqr_alpha


# 5. GRAFICAR LA DIVERSIDAD ALFA ---------------------------------------------

# Pasar la tabla a formato largo permite representar los cuatro índices
# en una misma figura. Cada panel conserva su propia escala vertical.

alpha_long <- reshape(
  alpha,
  varying = c(
    "observed",
    "shannon",
    "simpson_evenness",
    "chao1"
  ),
  v.names = "valor",
  timevar = "indice",
  times = c(
    "Riqueza observada",
    "Shannon",
    "Uniformidad de Simpson",
    "Chao1"
  ),
  direction = "long"
)


# Construir el gráfico.
# Las líneas grises unen las dos muestras del mismo animal.

grafico_alpha <- ggplot(
  alpha_long,
  aes(
    x = Sampling.position,
    y = valor
  )
) +
  geom_line(
    aes(group = animal),
    colour = "grey60",
    linewidth = 0.7
  ) +
  geom_boxplot(
    aes(colour = Sampling.position),
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Sampling.position),
    size = 3,
    position = position_jitter(
      width = 0.04,
      height = 0
    )
  ) +
  facet_wrap(
    ~ indice,
    scales = "free_y"
  ) +
  labs(
    title = "Diversidad alfa antes y después del chamuscado",
    x = "Etapa del procesamiento",
    y = "Valor del índice"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

grafico_alpha


# 6. DIVERSIDAD BETA ----------------------------------------------------------

# La diversidad beta compara la composición entre muestras.

# Bray-Curtis considera las abundancias relativas después de rarefacción.
dist_bray <- vegdist(
  t(assay(tse_d4, "rarefiedRelabundance")),
  method = "bray"
)

# Jaccard considera solamente presencia y ausencia de ASV.
dist_jaccard <- vegdist(
  t(assay(tse_d4, "pa")),
  method = "jaccard",
  binary = TRUE
)


# 7. PCoA --------------------------------------------------------------------

# PCoA representa las diferencias entre comunidades en dos dimensiones.
pcoa_bray <- ape::pcoa(dist_bray)
pcoa_jaccard <- ape::pcoa(dist_jaccard)

# Porcentaje de variación representado por los dos primeros ejes.
var_bray <- 100 * pcoa_bray$values$Relative_eig[1:2]
var_jaccard <- 100 * pcoa_jaccard$values$Relative_eig[1:2]

# Preparar las coordenadas de Bray-Curtis para graficar.
pcoa_bray_df <- data.frame(
  PCoA1 = pcoa_bray$vectors[, 1],
  PCoA2 = pcoa_bray$vectors[, 2],
  Sample = tse_d4$Sample,
  Sampling.position = tse_d4$Sampling.position
)

grafico_bray <- ggplot(
  pcoa_bray_df,
  aes(x = PCoA1, y = PCoA2, colour = Sampling.position)
) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -0.8, show.legend = FALSE) +
  labs(
    title = "PCoA basado en Bray-Curtis",
    x = paste0("PCoA1 (", round(var_bray[1], 1), " %)"),
    y = paste0("PCoA2 (", round(var_bray[2], 1), " %)"),
    colour = "Sitio"
  ) +
  theme_classic()

grafico_bray

# Preparar las coordenadas de Jaccard para graficar.
pcoa_jaccard_df <- data.frame(
  PCoA1 = pcoa_jaccard$vectors[, 1],
  PCoA2 = pcoa_jaccard$vectors[, 2],
  Sample = tse_d4$Sample,
  Sampling.position = tse_d4$Sampling.position
)

grafico_jaccard <- ggplot(
  pcoa_jaccard_df,
  aes(x = PCoA1, y = PCoA2, colour = Sampling.position)
) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -0.8, show.legend = FALSE) +
  labs(
    title = "PCoA basado en Jaccard",
    x = paste0("PCoA1 (", round(var_jaccard[1], 1), " %)"),
    y = paste0("PCoA2 (", round(var_jaccard[2], 1), " %)"),
    colour = "Sitio"
  ) +
  theme_classic()

grafico_jaccard


# 8. PERMANOVA ---------------------------------------------------------------

# PERMANOVA evalúa si la composición general difiere entre sitios.
set.seed(123)
permanova_bray <- adonis2(
  dist_bray ~ Sampling.position,
  data = as.data.frame(colData(tse_d4)),
  permutations = 999
)

set.seed(123)
permanova_jaccard <- adonis2(
  dist_jaccard ~ Sampling.position,
  data = as.data.frame(colData(tse_d4)),
  permutations = 999
)

permanova_bray
permanova_jaccard


# 9. PERMDISP ----------------------------------------------------------------

# PERMDISP comprueba si los grupos tienen una variabilidad interna similar.
disp_bray <- betadisper(dist_bray, tse_d4$Sampling.position)
disp_jaccard <- betadisper(dist_jaccard, tse_d4$Sampling.position)

set.seed(123)
permdisp_bray <- vegan::permutest(disp_bray, permutations = 999)

set.seed(123)
permdisp_jaccard <- vegan::permutest(disp_jaccard, permutations = 999)

permdisp_bray
permdisp_jaccard


# 10. NMDS (ACTIVIDAD OPCIONAL) -----------------------------------------------

# NMDS también representa las relaciones entre muestras, pero intenta conservar
# el orden de las disimilitudes en lugar de su magnitud exacta. Su calidad se
# evalúa mediante el stress: cuanto menor sea, mejor es la representación.
set.seed(123)
nmds_bray <- metaMDS(
  t(assay(tse_d4, "rarefiedRelabundance")),
  distance = "bray",
  k = 2,
  trymax = 100,
  autotransform = FALSE,
  trace = FALSE
)

nmds_bray$stress

# Extraer las coordenadas de las muestras.
nmds_coordenadas <- scores(nmds_bray, display = "sites")

nmds_df <- data.frame(
  NMDS1 = nmds_coordenadas[, 1],
  NMDS2 = nmds_coordenadas[, 2],
  Sample = tse_d4$Sample,
  Sampling.position = tse_d4$Sampling.position
)

grafico_nmds <- ggplot(
  nmds_df,
  aes(x = NMDS1, y = NMDS2, colour = Sampling.position)
) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -0.8, show.legend = FALSE) +
  labs(
    title = "NMDS basado en Bray-Curtis",
    subtitle = paste0("Stress = ", round(nmds_bray$stress, 3)),
    x = "NMDS1",
    y = "NMDS2",
    colour = "Sitio"
  ) +
  theme_classic()

grafico_nmds


# 11. GUARDAR RESULTADOS ------------------------------------------------------

dir.create("resultados_dia4", showWarnings = FALSE)

write.table(
  alpha,
  "resultados_dia4/diversidad_alfa.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  media_alpha,
  "resultados_dia4/medias_diversidad_alfa.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sd_alpha,
  "resultados_dia4/desvios_diversidad_alfa.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  resultados_wilcoxon,
  "resultados_dia4/tests_diversidad_alfa.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  "resultados_dia4/diversidad_alfa.png",
  grafico_alpha,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "resultados_dia4/PCoA_Bray_Curtis.png",
  grafico_bray,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "resultados_dia4/PCoA_Jaccard.png",
  grafico_jaccard,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "resultados_dia4/NMDS_Bray_Curtis.png",
  grafico_nmds,
  width = 7,
  height = 5,
  dpi = 300
)

saveRDS(tse_d4, "resultados_dia4/tse_dia4_final.rds")


# 12. PREGUNTAS PARA INTERPRETAR ---------------------------------------------

# 1. ¿Qué grupo presenta mayor riqueza observada?
# 2. ¿La diferencia de riqueza también aparece en Shannon y Simpson?
# 3. ¿Qué porcentaje de variación representan los dos primeros ejes del PCoA?
# 4. ¿Los sitios se separan visualmente en Bray-Curtis y Jaccard?
# 5. ¿Qué proporción de variación explica el sitio según PERMANOVA?
# 6. ¿El valor p de PERMANOVA indica evidencia de diferencias entre sitios?
# 7. ¿PERMDISP sugiere que los grupos tienen dispersiones diferentes?
# 8. ¿El NMDS muestra un patrón parecido al PCoA? ¿El stress es aceptable?
# 9. ¿Qué limitaciones tiene interpretar estos resultados con cuatro muestras
#    por grupo?
