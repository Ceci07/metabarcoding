#

## Load packages
library(mia)

## Beta diversity metrics like Bray-Curtis are often applied to relabundances
tse8_filt <- transformAssay(tse8_filt, assay.type="counts", method="relabundance")

# Other metrics like Aitchison to clr-transformed data
tse8_filt <- transformAssay(tse8_filt, assay.type="counts", method="clr", pseudocount=TRUE)

# Add group information Feces yes/no
tse$Group <- tse$SampleType == "Feces"


## Dissimilarity
## -------------
tse8_filt <- addDissimilarity(tse8_filt, assay.type = "counts" , method = "unifrac")


## Heatmap
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# Alternative
library(devtools)
install_github("jokergoo/ComplexHeatmap")

# load package
library(ComplexHeatmap)

# Annotation for samples
annotation <- HeatmapAnnotation(sample_type = tse8_filt[["Sampling.position"]])

# Create a heatmap
Heatmap(
  metadata(tse8_filt)[["unifrac"]],
  heatmap_legend_param = list(title = "Unifrac"),
  bottom_annotation = annotation
)


## Divergence
## ----------
tse8_filt <- addDivergence(tse8_filt, assay.type = "counts", reference = "median", FUN = getDissimilarity)


## PERMANOVA
## ---------
res <- getPERMANOVA(
    tse8_filt,
    assay.type = "relabundance",
    formula = x ~ Sampling.position
    )
res


## Unsupervised ordination
## -----------------------
library(scater)

install.packages("ggExtra")
library(ggExtra)


# Run PCoA on relabundance assay with Bray-Curtis distances
tse8_filt <- addMDS(
    tse8_filt,
    FUN = getDissimilarity,
    method = "bray",
    assay.type = "relabundance",
    name = "MDS_bray")

## Create ggplot object
p <- plotReducedDim(tse8_filt, "MDS_bray", colour_by = "Sampling.position")

# Calculate explained variance
e <- attr(reducedDim(tse8_filt, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p + labs(
    x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
    y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")
    )
p

#
p <- ggMarginal(p, type = "boxplot", groupFill = TRUE)
p

# Run NMDS on relabundance assay with Bray-Curtis distances
tse8_filt <- addNMDS(
    tse8_filt,
    FUN = getDissimilarity,
    method = "bray",
    assay.type = "relabundance",
    name = "NMDS_bray")

## Run MDS on clr assay with Aitchison distances
tse8_filt <- addMDS(
    tse8_filt,
    FUN = getDissimilarity,
    method = "euclidean",
    assay.type = "clr",
    name = "MDS_aitchison")

## Run NMDS on clr assay with Euclidean distances
tse8_filt <- addNMDS(
    tse8_filt,
    FUN = getDissimilarity,
    method = "euclidean",
    assay.type = "clr",
    name = "NMDS_aitchison")


## Load package for multi-panel plotting
library(patchwork)

# Generate plots for all 4 reducedDims
plots <- lapply(
    c("MDS_bray", "MDS_aitchison", "NMDS_bray", "NMDS_aitchison"),
    plotReducedDim,
    object = tse8_filt,
    colour_by = "Sampling.position")

# Generate multi-panel plot
wrap_plots(plots) +
    plot_layout(guides = "collect")


## Unifrac
tse8_filt <- addMDS(
    tse8_filt,
    FUN = getDissimilarity,
    method = "unifrac",
    name = "unifrac",
    tree = rowTree(tse8_filt),
    ntop = nrow(tse8_filt),
    assay.type = "counts")

plotReducedDim(tse8_filt, "unifrac", colour_by = "Sampling.position")

## Rarefaction to mitigate impacts of uneven sequencing effort
# Calculate the list of sequencing depths across samples
sequencing_depths <- colSums(assay(tse8_filt))

# Calculate variation between highest and lowest sequencing depth
depth_variation <- max(sequencing_depths)/min(sequencing_depths)
depth_variation
[1] 2.258765


## Centered log-ratio transformation to properly apply Aitchison distance
tse8_filt <- transformAssay(
    tse8_filt,
    assay.type = "counts",
    method = "clr",
    pseudocount = 1
    )

## Run MDS on clr assay with Aitchison distance
tse8_filt <- addMDS(
    tse8_filt,
    FUN = getDissimilarity,
    method = "euclidean",
    assay.type = "clr",
    name = "MDS_aitchison"
    )

## with rarefaction
# Define custom transformation function..
clr <- function (x) {
    vegan::decostand(x, method = "clr", pseudocount = 1)
}

# Run transformation after rarefactions before calculating the beta diversity
tse8_filt <- addMDS(
    tse8_filt,
    assay.type = "counts",
    FUN = getDissimilarity,
    method = "euclidean",
    niter = 10, # Number of iterations
    sample = min(colSums(assay(tse8_filt, "counts"))), # Rarefaction depth
    transf = clr, # Applied transformation
    replace = TRUE,
    name = "MDS_aitchison_rarefied"
    )

## Generate plots for non-rarefied and rarefied Bray-Curtis distances scaled by
plots <- lapply(
    c("MDS_aitchison", "MDS_aitchison_rarefied"),
    plotReducedDim,
    object = tse8_filt,
    colour_by = "Sampling.position"
    )

# Generate multi-panel plot
wrap_plots(plots) +
    plot_layout(guides = "collect")

## Plot the correlation between principal coordinates for both the rarefied and non-rarefied distance
install.packages("ggpubr")
library(ggpubr)

p <- lapply(1:2, function(i){
    # Principal axes are sign invariant;
    # align the signs; if the correlation is negative then swap the other axis
    original <- reducedDim(tse8_filt, "MDS_aitchison")[, i]
    rarefied <- reducedDim(tse8_filt, "MDS_aitchison_rarefied")[, i]
  
    temp  <- ggplot(
        data = data.frame(original, rarefied),
        aes(x = original, y = rarefied)
        ) +
        geom_point() + geom_smooth(method = "lm") +
        stat_cor(method = "pearson") +
        labs(title = paste0("Principal coordinate ", i))
    return(temp)
})
wrap_plots(p)


## Differential abundance
## ----------------------
# ANCOM-BC

# Load package
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
BiocManager::install("microbiome")

install.packages("knitr")
install.packages("dplyr")

library(ANCOMBC) 
library(microbiome)
library(knitr)
library(dplyr)

## Run ANCOM-BC at the genus level and only including the prevalent genera
ancombc2_out <- ancombc2(
    data = tse8_filt,
    assay.type = "counts",
    fix_formula = "Sampling.position",
    p_adj_method = "fdr",
    prv_cut = 0,
    group = "Sampling.position",
    struc_zero = TRUE,
    neg_lb = TRUE,
    # multi group comparison is deactivated automatically
    global = TRUE)

# store the FDR adjusted results
ancombc2_out$res |>
  dplyr::select(taxon, lfc_Sampling.positionSticking, q_Sampling.positionSticking) |>
  filter(q_Sampling.positionSticking < 0.05) |>
  arrange(q_Sampling.positionSticking) |>
  head() |>
  kable()

|taxon                       | lfc_Sampling.positionSticking| q_Sampling.positionSticking|
|:---------------------------|-----------------------------:|---------------------------:|
|Moraxella                   |                      3.034162|                   0.0283150|
|Chryseobacterium            |                      3.793471|                   0.0310717|
|Bergeyella                  |                      3.373888|                   0.0310717|
|Flavobacterium              |                      4.333129|                   0.0310717|
|Clostridium sensu stricto 1 |                      1.882225|                   0.0310717|
|Corynebacterium             |                      3.299914|                   0.0375974|

