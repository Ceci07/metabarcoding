# ============================================================
# CURSO DE METABARCODING - DÍA 3
# mia + TreeSummarizedExperiment + exploración y QC
#
# Dataset: Zwirzitz et al. (2020)
# The sources and transmission routes of microbial populations
# throughout a meat processing facility
# npj Biofilms and Microbiomes 6:26
# DOI: 10.1038/s41522-020-0136-z
# ENA: PRJEB37434
#
# En este curso usamos SOLO el conjunto Illumina.
#
# OBJETIVOS
# 1. Entender la estructura de TreeSummarizedExperiment (TreeSE)
# 2. Integrar/explorar metadata
# 3. Explorar profundidad de secuenciación
# 4. Explorar dominancia y prevalencia
# 5. Discutir qué muestras/features tiene sentido filtrar
# 6. Introducir aglomeración, transformación y composición
#
# INPUTS ESPERADOS
# data/sampleData_alimentos.csv
# data/tse_illumina.rds
#
# El RDS debe contener:
# - assay(tse, "counts") : matriz ASV x muestra
# - rowData(tse)        : taxonomía
# - rowTree(tse)        : árbol (si está disponible)
# - colData(tse)        : metadata, o bien incorporar el CSV
# ============================================================

library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(tidyverse)

# ============================================================
# 1. METADATA
# ============================================================

metadata <- read.csv(
    "data/sampleData_alimentos.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

head(metadata)
str(metadata)

# Variables:
# Run, Sample, Sampling position, Farmer

table(metadata[["Sampling position"]])

# ============================================================
# 2. CARGAR TREESE
# ============================================================

# Este objeto es el resultado del procesamiento realizado
# en el Día 2 (DADA2 + taxonomía + árbol).

library(phyloseq)

ps<-readRDS("Downloads/phyloseq.RDS")

# Convertimos phyloseq a TreeSummarizedExperiment (TSE) via mia
tse <- convertFromPhyloseq(ps)


class(tse)
tse
dim(tse)

# ============================================================
#  Chequeo de algunas secuencias ruido
# ============================================================

euka<-which(rowData(tse)$Kingdom=="Eukaryota")
length(euka)

#tse<-tse[-euka,]
#dim(tse)


mito<-which(rowData(tse)$Family=="Mitochondria")
length(mito)

tse<-tse[-mito,]
dim(tse)

cloro<-which(rowData(tse)$Family=="Chloroplast")
length(cloro)

#tse<-tse[-cloro,]
#dim(tse)

# ============================================================
# 3. ANATOMÍA DEL TREESE
# ============================================================

# Samples = columnas
# Features/ASVs = filas

nrow(tse)
ncol(tse)

head(colnames(tse))
head(rownames(tse))

# ASSAYS: abundancias
assayNames(tse)
assays(tse)

counts <- assay(tse, "counts")
dim(counts)
counts[1:5, 1:5]

# COLDATA: metadata de las muestras
colData(tse)

# ROWDATA: metadata de features / taxonomía
rowData(tse)
colnames(rowData(tse))

# Árbol filogenético
rowTree(tse)
rowLinks(tse)

# ============================================================
# 4. INCORPORAR / COMPROBAR METADATA
# ============================================================

# Los identificadores deben coincidir.
head(colnames(tse))
head(metadata$Run)

metadata2 <- metadata[match(colnames(tse), metadata$Run), ]

# TRUE = el orden coincide correctamente
all(metadata2$Run == colnames(tse))

# Si el TreeSE no contiene estas variables, incorporarlas:
colData(tse)$Run <- metadata2$Run
colData(tse)$Sample <- metadata2$Sample
colData(tse)$`Sampling position` <- metadata2$`Sampling position`
colData(tse)$Farmer <- metadata2$Farmer

# Explorar las categorías
table(tse$`Sampling position`)

# ============================================================
# 5. TAXONOMÍA
# ============================================================

taxonomyRanks(tse)
checkTaxonomy(tse)

# Ejemplo: visualizar algunos features
rowData(tse)[1:10, taxonomyRanks(tse)]

# ============================================================
# 6. SUBSETTING: UNA PROPIEDAD FUNDAMENTAL
# ============================================================

# Al hacer subset, las piezas del objeto permanecen coordinadas.

tse_small <- tse[, 1:5]

dim(tse_small)
colData(tse_small)
rowData(tse_small)[1:5, ]

# ============================================================
# 7. PROFUNDIDAD DE SECUENCIACIÓN
# ============================================================

library_size <- colSums(assay(tse, "counts"))

summary(library_size)

tse$library_size <- library_size

ggplot(
    as.data.frame(colData(tse)),
    aes(x = library_size)
) +
    geom_histogram(bins = 30) +
    theme_bw() +
    labs(
        x = "Número de reads",
        y = "Número de muestras",
        title = "Profundidad de secuenciación"
    )

# Por tipo de muestra
ggplot(
    as.data.frame(colData(tse)),
    aes(
        x = `Sampling position`,
        y = library_size
    )
) +
    geom_boxplot() +
    geom_jitter(width = 0.15, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        x = "Sampling position",
        y = "Número de reads"
    )

# PREGUNTAS:
# - ¿Hay muestras con muy poca profundidad?
# - ¿Hay diferencias de profundidad entre grupos?
# - ¿Una muestra con poca profundidad es necesariamente mala?
#
# NO necesariamente. La profundidad debe interpretarse
# junto con la pregunta y el resto del QC.




# ============================================================
# 8. DOMINANCIA
# ============================================================

# Dominancia pregunta cuánto contribuyen los features
# a las lecturas de una muestra.

dominance <- summarizeDominance(
    tse,
    assay.type = "counts"
)

head(dominance)

# También podemos preparar abundancias relativas.
tse_rel <- transformAssay(
    tse,
    assay.type = "counts",
    method = "relabundance"
)

# ============================================================
# 9. PREVALENCIA
# ============================================================

# Prevalencia = en cuántas muestras aparece un feature.

tse <- addPrevalence(
    tse,
    assay.type = "counts",
    detection = 0
)

head(rowData(tse)$prevalence)

ggplot(
    as.data.frame(rowData(tse)),
    aes(x = prevalence)
) +
    geom_histogram(bins = 30) +
    theme_bw() +
    labs(
        x = "Número de muestras en las que aparece el feature",
        y = "Número de features"
    )

# ============================================================
# 10. Curvas de rarefacción, riqueza de ASVs y filtrado por prevalencia
# ============================================================

# Extraemos la matriz de conteos para rarefaction curves
count_matrix <- t(assay(tse, "counts"))

library(vegan)
rarecurve(count_matrix, step = 100, label = TRUE, cex = 0.5)

# Conteo de ASVs por muestra y estadísticos resumen
asv_per_sample <- apply(counts(tse), 2, function(x) sum(x > 0))
asv_summary <- c(min = min(asv_per_sample), mean = mean(asv_per_sample), max = max(asv_per_sample))

print(asv_per_sample)

print(asv_summary)


# conservar features presentes en al menos 10% de las muestras.

tse_prev10 <- subsetByPrevalent(
    tse,
    prevalence = 0.10,
    detection = 0,
    update.tree = TRUE
)

dim(tse)
dim(tse_prev10)

# ¿Cuántos features se eliminaron?
nrow(tse) - nrow(tse_prev10)

# IMPORTANTE:
# 10% NO es un umbral universal.
# Es sólo un ejemplo para enseñar el procedimiento.
#
# La decisión depende de:
# - objetivo del estudio
# - tamaño de muestra
# - interés en taxones raros
# - ruido esperado
# - diseño experimental
#
# En un análisis real conviene evaluar la robustez
# de las conclusiones frente a distintos filtros.

# ============================================================
# 11. ¿FILTRAR MUESTRAS POR PROFUNDIDAD?
# ============================================================

summary(tse$library_size)

# NO eliminar automáticamente la muestra de menor profundidad.
#
# Si se decide usar un umbral, debe explicitarse:
#
# min_reads <- 5000
# tse_depth <- tse[, tse$library_size >= min_reads]
#
# Después comparar:
# dim(tse)
# dim(tse_depth)
#
# Pregunta:
# ¿la muestra excluida es un problema técnico
# o una observación biológicamente importante?


##############################################################################
# Rarefacción
##############################################################################


min_sample_sums<-min(apply(counts(tse), 2,sum))
tse <- rarefyAssay(tse, assay.type = "counts", sample = min_sample_sums, replace = FALSE)

# Obtener abundancias relativas (útil para comparaciones y visualizaciones)
tse <- transformAssay(tse, assay.type = "subsampled", method = "relabundance", name="rarefiedRelabundance")

tse <- transformAssay(tse, method = "pa")

# Conteo de ASVs por muestra y estadísticos resumen
asv_per_sample <- apply(assay(tse, "subsampled"), 2, function(x) sum(x > 0))
asv_summary <- c(min = min(asv_per_sample), mean = mean(asv_per_sample), max = max(asv_per_sample))

print(asv_per_sample)
   
   
print(asv_summary)


# Exportar tablas crudas y relativas si se desean compartir en repositorio
write.table(assay(tse, "counts"), file = "conteos_ASVs_prevalencia050.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(assay(tse, "relabundance"), file = "abundanciaRel_ASVs_prevalencia05.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(assay(tse, "pa"), file = "pa_ASVs_prevalencia05.txt", sep = "\t", quote = FALSE, col.names = NA)
write.table(rowData(tse), file = "taxonomiaASVs.txt", sep = "\t", quote = FALSE, col.names = NA)



# ============================================================
# 12. AGLOMERACIÓN TAXONÓMICA
# ============================================================

# ASV = máxima resolución de la secuencia observada.
# Para algunas visualizaciones podemos agrupar por género.

tse_genus <- agglomerateByRank(
    tse_prev10,
    rank = "Genus"
)

dim(tse_genus)

# Aglomerar NO es lo mismo que filtrar:
# - filtrado = eliminar features/muestras
# - aglomeración = agrupar features

# ============================================================
# 13. TRANSFORMACIÓN
# ============================================================

# Abundancia relativa:
tse_genus_rel <- transformAssay(
    tse_genus,
    assay.type = "counts",
    method = "relabundance"
)

assayNames(tse_genus_rel)

# Para algunos análisis exploratorios puede usarse CLR.
# Es una representación derivada, no reemplaza los counts originales.

tse_clr <- transformAssay(
    tse_genus,
    assay.type = "counts",
    method = "clr",
    pseudocount = TRUE
)


##############################################################################
# Porcentaje de lecturas clasificadas a distintos niveles taxonómicos
##############################################################################

# Calculamos la fracción de abundancia que quedó sin asignar en cada rango taxonómico
# 1 - mean(colSums(relabundance de los features NA en ese nivel))  => fracción asignada
prop_assigned<-vector()
prop_assigned <- rbind(
  Kingdom = 1 - apply(assay(tse[is.na(rowData(tse)$Kingdom),], "rarefiedRelabundance"),2,sum),
  Phylum  = 1 - apply(assay(tse[is.na(rowData(tse)$Phylum),], "rarefiedRelabundance"),2,sum),
  Class   = 1 - apply(assay(tse[is.na(rowData(tse)$Class),], "rarefiedRelabundance"),2,sum),
  Order   = 1 - apply(assay(tse[is.na(rowData(tse)$Order),], "rarefiedRelabundance"),2,sum),
  Family  = 1 - apply(assay(tse[is.na(rowData(tse)$Family),], "rarefiedRelabundance"),2,sum),
  Genus   = 1 - apply(assay(tse[is.na(rowData(tse)$Genus),], "rarefiedRelabundance"),2,sum),
  Species = 1 - apply(assay(tse[is.na(rowData(tse)$Species),], "rarefiedRelabundance"),2,sum)
)
print(prop_assigned)


# Número de categorias observadas en cada nivel
n_taxa_levels <- c(
  Kingdom = length(unique(rowData(tse)$Kingdom)),
  Phylum  = length(unique(rowData(tse)$Phylum)),
  Class   = length(unique(rowData(tse)$Class)),
  Order   = length(unique(rowData(tse)$Order)),
  Family  = length(unique(rowData(tse)$Family)),
  Genus   = length(unique(rowData(tse)$Genus))
)

print(n_taxa_levels)


####################################################################

library(ComplexHeatmap)

# Annotation for samples
annotation <- HeatmapAnnotation(sample_type = colData(tse)$"Sampling.position")

# Create a heatmap

pdf("heatmap_jaccard.pdf")
	Heatmap(
    		metadata(tse)[["jaccard"]],
    		heatmap_legend_param = list(title = "jaccard"),
    		bottom_annotation = annotation
	)
dev.off()


pdf("heatmap_jaccard.pdf")
	Heatmap(
    		metadata(tse)[["jaccard"]],
    		heatmap_legend_param = list(title = "jaccard"),
    		bottom_annotation = annotation
	)
dev.off()

##############################################################################
#  Community composition
##############################################################################

tse_phylum<-agglomerateByRank(tse, rank="Phylum", empty.rm=F, na.rm=F)


top_phyla <- getTop(tse_phylum, top= 10, assay_name="counts",na.rm=F)
top_phyla<-rowData(tse_phylum[top_phyla,])$Phylum

# Reetiquetar
tax <- rowData(tse_phylum)

tax$Phylum2 <- ifelse(
  tax$Phylum %in% top_phyla,
  as.character(tax$Phylum),
  "Other"
)

rowData(tse_phylum)$Phylum2 <- tax$Phylum2

tse_phylum<- transformAssay(tse_phylum, method = "relabundance")


p4 <- plotAbundance(
    tse_phylum, group="Phylum2", order.row.by="abund", assay.type="relabundance",
    assay = "relabundance", order.col.by="km") + theme(axis.text.x = element_text(angle = 90))

pdf("filosMasAbundantes.pdf")
	p4
dev.off()





# ============================================================
# 14. COMPOSICIÓN: EJEMPLO CON NUESTROS DATOS
# ============================================================

# Elegimos tres tipos de muestra que permiten una discusión biológica:
# - Sticking: piel/animal
# - Anal swab: hisopado anal
# - Classification gloves: guantes de clasificación

keep_positions <- c(
    "Sticking",
    "Anal swab",
    "Classification gloves"
)

tse_demo <- tse_genus_rel[
    ,
    tse_genus_rel$`Sampling position` %in% keep_positions
]

# Visualizar abundancias relativas
plotAbundance(
    tse_demo,
    assay.type = "relabundance",
    group = "Sampling position"
)

# PREGUNTAS PARA DISCUSIÓN:
#
# 1. ¿Esperamos que piel y guantes sean similares?
# 2. ¿Esperamos que el hisopado anal sea diferente?
# 3. ¿Qué taxones parecen caracterizar cada tipo?
# 4. ¿Una diferencia visual implica necesariamente
#    una diferencia estadísticamente significativa?
# 5. ¿Podría la profundidad de secuenciación explicar
#    alguna diferencia?



# ============================================================
# 16. CHECKLIST DEL DÍA 3
# ============================================================

# Debemos poder responder:

# ¿Cuántas muestras?
dim(tse)

# ¿Cuántos features?
nrow(tse)

# ¿Qué metadata tenemos?
colnames(colData(tse))

# ¿Qué taxonomía tenemos?
taxonomyRanks(tse)

# ¿Tenemos árbol?
rowTree(tse)

# ¿Cómo es la profundidad?
summary(tse$library_size)

# ¿Cómo es la prevalencia?
summary(rowData(tse)$prevalence)

# ¿Qué criterio de filtrado elegimos?
#
# ¿Qué información perdemos al filtrar?
#
# ¿Qué pregunta biológica queremos responder?

# ============================================================
# 17. GUARDAR EL OBJETO PARA EL DÍA 4
# ============================================================

# Cuando el grupo haya decidido el filtro final:
#
# tse_D3 <- tse_prev10
#
# saveRDS(
#     tse_D3,
#     "data/tse_illumina_D3_filtered.rds"
# )

# ============================================================
# DÍA 4
# ============================================================
#
# El objeto final de hoy será la entrada para:
# - alpha diversity
# - beta diversity
# - distancias
# - ordenación
# - visualización
#
# ============================================================
