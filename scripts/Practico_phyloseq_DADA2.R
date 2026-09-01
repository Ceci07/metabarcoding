# ==============================================================================
# PRACTICO DE PHYLOSEQ A PARTIR DE DADA2
# ==============================================================================

library(phyloseq)
library(ggplot2)
library(dplyr)

set.seed(123)

# ------------------------------------------------------------------------------
# 1. ARCHIVOS DE ENTRADA
# ------------------------------------------------------------------------------

data_dir <- "."

seqtab_file <- file.path(data_dir, "seqtab.nochim.rds")
taxa_file   <- file.path(data_dir, "taxa_gtdb.gz")
meta_file   <- file.path(data_dir, "sampleData_alimentos.csv")

# ------------------------------------------------------------------------------
# 2. TABLA DE ASVs DE DADA2
# ------------------------------------------------------------------------------

seqtab.nochim <- readRDS(seqtab_file)

class(seqtab.nochim)
dim(seqtab.nochim)

head(rownames(seqtab.nochim))
substr(colnames(seqtab.nochim)[1:3], 1, 80)

# PREGUNTA:
# ¿Que representan las filas y las columnas de seqtab.nochim?
#
# PREGUNTA:
# ¿Cual es la diferencia entre una OTU y una ASV?

# ------------------------------------------------------------------------------
# 3. METADATOS
# ------------------------------------------------------------------------------

metadata <- read.csv(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

head(metadata)
dim(metadata)

# PREGUNTA:
# ¿Que columnas tiene el archivo de metadatos?
colnames(metadata)

# PREGUNTA:
# ¿Que columna puede utilizarse como identificador de las muestras?

rownames(metadata) <- metadata$Run

# Creamos un nombre sin espacios para usarlo mas facilmente en el practico.
metadata$Sampling_position <- metadata[["Sampling position"]]

# ------------------------------------------------------------------------------
# 4. TAXONOMIA
# ------------------------------------------------------------------------------

load(taxa_file)

class(taxa_gtdb)
dim(taxa_gtdb)
head(taxa_gtdb)

# PREGUNTA:
# ¿Que niveles taxonomicos tiene la tabla de taxas?
colnames(taxa_gtdb)

# ------------------------------------------------------------------------------
# 5. COMPROBAR LA COINCIDENCIA ENTRE OBJETOS
# ------------------------------------------------------------------------------

# Muestras
samples_in_counts <- rownames(seqtab.nochim)
samples_in_meta   <- rownames(metadata)

setdiff(samples_in_counts, samples_in_meta)
setdiff(samples_in_meta, samples_in_counts)

# Reordenar metadata de acuerdo con la tabla de ASVs
metadata_ps <- metadata[samples_in_counts, , drop = FALSE]

identical(
  rownames(seqtab.nochim),
  rownames(metadata_ps)
)

# ASVs
asvs_in_counts <- colnames(seqtab.nochim)
asvs_in_tax    <- rownames(taxa_gtdb)

length(asvs_in_counts)
length(asvs_in_tax)

common_asvs <- intersect(
  asvs_in_counts,
  asvs_in_tax
)

length(common_asvs)

# Esta comprobacion no es estrictamente necesaria si confiamos en que las ASVs
# y las muestras coinciden entre los objetos. Sin embargo, la falta de
# coincidencia de IDs de muestras o ASVs es una de las fuentes mas comunes de
# error al construir objetos phyloseq.

seqtab_ps <- seqtab.nochim[, common_asvs, drop = FALSE]
taxa_ps   <- taxa_gtdb[common_asvs, , drop = FALSE]

# ------------------------------------------------------------------------------
# 6. RENOMBRAR LAS SECUENCIAS COMO ASV1, ASV2, ..., ASVn
# ------------------------------------------------------------------------------

# Guardamos primero la correspondencia entre el nuevo ID y la secuencia original.
asv_sequences <- colnames(seqtab_ps)
asv_ids <- paste0("ASV", seq_along(asv_sequences))

asv_map <- data.frame(
  ASV = asv_ids,
  Sequence = asv_sequences,
  stringsAsFactors = FALSE
)

# Renombramos de la misma forma la tabla de abundancias y la tabla taxonomica.
colnames(seqtab_ps) <- asv_ids
rownames(taxa_ps) <- asv_ids

head(asv_map)
head(colnames(seqtab_ps))
head(rownames(taxa_ps))

# ------------------------------------------------------------------------------
# 7. CREAR LOS COMPONENTES DE PHYLOSEQ
# ------------------------------------------------------------------------------

# Tabla de abundancias
OTU <- otu_table(
  seqtab_ps,
  taxa_are_rows = FALSE
)

OTU
taxa_are_rows(OTU)

# Tabla taxonomica
TAX <- tax_table(
  as.matrix(taxa_ps)
)

TAX

# Metadatos
SAM <- sample_data(
  metadata_ps
)

SAM

# ------------------------------------------------------------------------------
# 8. CREAR EL OBJETO PHYLOSEQ
# ------------------------------------------------------------------------------

ps <- phyloseq(
  otu_table = OTU,
  tax_table = TAX,
  sample_data = SAM
)

ps

# PREGUNTA:
# ¿Cuantas muestras y cuantas ASVs contiene el objeto phyloseq?

nsamples(ps)
ntaxa(ps)

# ------------------------------------------------------------------------------
# 9. EXPLORAR EL OBJETO
# ------------------------------------------------------------------------------

head(sample_names(ps))
head(taxa_names(ps))

rank_names(ps)
sample_variables(ps)

otu_table(ps)[1:5, 1:5]
head(as(tax_table(ps), "matrix"))
head(data.frame(sample_data(ps)))

# ------------------------------------------------------------------------------
# 10. PROFUNDIDAD DE SECUENCIACION
# ------------------------------------------------------------------------------

depth <- sample_sums(ps)

head(depth)
summary(depth)

min(depth)
median(depth)
max(depth)

# PREGUNTA:
# ¿Cuantas veces mas reads tiene la muestra con mayor profundidad respecto de
# la muestra con menor profundidad?
max(depth) / min(depth)

depth_df <- data.frame(
  Reads = sample_sums(ps),
  sample_data(ps),
  check.names = FALSE
)

ggplot(
  depth_df,
  aes(
    x = reorder(Sample, Reads),
    y = Reads,
    fill = Sampling_position
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Muestra",
    y = "Numero de reads",
    fill = "Posicion de muestreo",
    title = "Profundidad de secuenciacion"
  ) +
  theme_bw()

# Abundancia total de cada ASV
head(
  sort(
    taxa_sums(ps),
    decreasing = TRUE
  ),
  20
)

# ------------------------------------------------------------------------------
# 11. SUBSETTING Y FILTRADO
# ------------------------------------------------------------------------------

# Filtrar muestras
ps_sticking <- subset_samples(
  ps,
  Sampling_position == "Sticking"
)

ps_sticking

# Filtrar un grupo taxonomico
ps_firmicutes <- subset_taxa(
  ps,
  Phylum == "Firmicutes"
)

ps_firmicutes

# Mantener ASVs con al menos 100 reads en todo el dataset
ps_abundant <- prune_taxa(
  taxa_sums(ps) >= 100,
  ps
)

ps_abundant

# Mantener ASVs presentes en al menos 5 muestras
ps_prevalent <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= 5,
  prune = TRUE
)

ps_prevalent

# PREGUNTA:
# ¿Cual es la diferencia entre abundancia y prevalencia?
# ¿Puede un taxon ser muy abundante pero tener baja prevalencia?
# ¿Puede tener baja abundancia pero alta prevalencia?

# ------------------------------------------------------------------------------
# 12. AGLOMERACION TAXONOMICA CON tax_glom()
# ------------------------------------------------------------------------------

ps_phylum <- tax_glom(
  ps,
  taxrank = "Phylum",
  NArm = TRUE
)

ps_genus <- tax_glom(
  ps,
  taxrank = "Genus",
  NArm = TRUE
)

c(
  ASVs = ntaxa(ps),
  Phyla = ntaxa(ps_phylum),
  Genera = ntaxa(ps_genus)
)

# PREGUNTA:
# ¿Que hace tax_glom() con las abundancias de ASVs que pertenecen al mismo
# grupo taxonomico?

# ------------------------------------------------------------------------------
# 13. COMPOSICION TAXONOMICA
# ------------------------------------------------------------------------------

ps_phylum_rel <- transform_sample_counts(
  ps_phylum,
  function(x) x / sum(x)
)

plot_bar(
  ps_phylum_rel,
  x = "Sample",
  fill = "Phylum"
) +
  labs(
    x = "Muestra",
    y = "Abundancia relativa",
    title = "Composicion taxonomica a nivel de Phylum"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

# ------------------------------------------------------------------------------
# 14. ¿POR QUE NORMALIZAR?
# ------------------------------------------------------------------------------

# Ejemplo:
#
# Muestra T1 = 20,000 reads
# Muestra T2 = 200,000 reads
#
# T2 logro secuenciar 10 veces mas reads y, por lo tanto, contiene 10 veces mas
# observaciones de taxa. Esta diferencia puede deberse al proceso tecnico de
# secuenciacion y no implica que T2 contenga 10 veces mas microorganismos.
#
# Por eso los counts brutos entre muestras con distinta profundidad no deben
# interpretarse directamente como abundancias microbianas absolutas.

# ------------------------------------------------------------------------------
# 15. DATOS COMPOSICIONALES
# ------------------------------------------------------------------------------

# Ejemplo conceptual:
#
#                 T1      T2
# Taxon A         100     200
# Taxon B         100     100
#
# En abundancia absoluta, B no cambio.
#
# Pero al convertir a proporciones:
#
#                 T1       T2
# Taxon A         50%      66.7%
# Taxon B         50%      33.3%
#
# La abundancia relativa de B disminuye aunque su abundancia absoluta no haya
# cambiado.
#
# PREGUNTA:
# ¿Por que ocurre esta disminucion aparente del Taxon B?

# ------------------------------------------------------------------------------
# 16. ABUNDANCIA RELATIVA
# ------------------------------------------------------------------------------

ps_rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

head(sample_sums(ps_rel))

# PREGUNTA:
# ¿Cuanto suma cada muestra despues de transformar a abundancia relativa?

# ------------------------------------------------------------------------------
# 17. CPM: COUNTS PER MILLION
# ------------------------------------------------------------------------------

ps_cpm <- transform_sample_counts(
  ps,
  function(x) (x / sum(x)) * 1e6
)

head(sample_sums(ps_cpm))

# PREGUNTA:
# ¿Cuanto suma cada muestra en CPM?
# ¿Que porcentaje representa un taxon con 50,000 CPM?

# ------------------------------------------------------------------------------
# 18. RAREFACCION
# ------------------------------------------------------------------------------

min_depth <- min(sample_sums(ps))
min_depth

ps_rare <- rarefy_even_depth(
  ps,
  sample.size = min_depth,
  rngseed = 123,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

sample_sums(ps_rare)

depth_rare_df <- data.frame(
  Reads = sample_sums(ps_rare),
  sample_data(ps_rare),
  check.names = FALSE
)

ggplot(
  depth_rare_df,
  aes(
    x = Sample,
    y = Reads,
    fill = Sampling_position
  )
) +
  geom_col() +
  labs(
    x = "Muestra",
    y = "Numero de reads",
    fill = "Posicion de muestreo",
    title = "Profundidad despues de rarefaccion"
  ) +
  theme_bw()

# PREGUNTA:
# ¿Que diferencia observamos entre la profundidad de las muestras antes y
# despues de rarefactar?

# ------------------------------------------------------------------------------
# 19. DIVERSIDAD ALFA
# ------------------------------------------------------------------------------

# Observed = numero de ASVs observados
# Shannon = considera riqueza y equidad
# Simpson = da mayor peso relativo a taxa dominantes

alpha <- estimate_richness(
  ps_rare,
  measures = c(
    "Observed",
    "Shannon",
    "Simpson"
  )
)

alpha$Run <- rownames(alpha)

meta_rare <- data.frame(sample_data(ps_rare))
meta_rare$Run <- rownames(meta_rare)

alpha_df <- merge(
  alpha,
  meta_rare,
  by = "Run",
  all.x = TRUE
)

head(alpha_df)

# Pasamos los tres indices a formato largo usando solamente R base.
alpha_long <- data.frame(
  Sampling_position = rep(alpha_df$Sampling_position, times = 3),
  Metric = factor(
    rep(
      c("Observed", "Shannon", "Simpson"),
      each = nrow(alpha_df)
    ),
    levels = c("Observed", "Shannon", "Simpson")
  ),
  Value = c(
    alpha_df$Observed,
    alpha_df$Shannon,
    alpha_df$Simpson
  )
)

ggplot(
  alpha_long,
  aes(
    x = Sampling_position,
    y = Value,
    fill = Sampling_position
  )
) +
  geom_boxplot() +
  facet_wrap(
    ~Metric,
    scales = "free_y"
  ) +
  labs(
    x = "Posicion de muestreo",
    y = "Indice de diversidad"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# PREGUNTAS:
# ¿Que diferencia conceptual hay entre riqueza, Shannon y Simpson?
# ¿Dos muestras con la misma riqueza pueden tener distinto Shannon?

# ------------------------------------------------------------------------------
# 20. DIVERSIDAD BETA Y PCoA
# ------------------------------------------------------------------------------

# PCoA (Principal Coordinates Analysis) parte de una matriz de distancias entre
# muestras y busca representarlas en pocos ejes. Cada punto representa una
# muestra; dos puntos proximos corresponden a comunidades mas similares segun
# la metrica de distancia utilizada.

# Bray-Curtis: considera diferencias de abundancia.
dist_bray <- distance(
  ps_rel,
  method = "bray"
)

ord_bray <- ordinate(
  ps_rel,
  method = "PCoA",
  distance = dist_bray
)

plot_ordination(
  ps_rel,
  ord_bray,
  color = "Sampling_position"
) +
  geom_point(size = 3) +
  labs(
    color = "Posicion de muestreo",
    title = "PCoA - Bray-Curtis"
  ) +
  theme_bw()

# Jaccard binario: considera presencia/ausencia y en base a eso computa distancias.
dist_jaccard <- distance(
  ps_rare,
  method = "jaccard",
  binary = TRUE
)

# La PCoA utiliza la matriz de distancias de Jaccard para obtener coordenadas
# que representen las relaciones entre las muestras en pocos ejes.
ord_jaccard <- ordinate(
  ps_rare,
  method = "PCoA",
  distance = dist_jaccard
)

# Plot de las coordenadas de cada muestra en los primeros ejes de la PCoA.
plot_ordination(
  ps_rare,
  ord_jaccard,
  color = "Sampling_position"
) +
  geom_point(size = 3) +
  labs(
    color = "Posicion de muestreo",
    title = "PCoA - Jaccard"
  ) +
  theme_bw()

# PREGUNTAS:
# ¿Que representa cada punto de una PCoA?
# ¿Que significa que dos puntos esten proximos?
# ¿Por que Bray-Curtis y Jaccard pueden producir ordenaciones diferentes?
# ¿Una separacion visual en una PCoA demuestra una diferencia estadisticamente
# significativa?
