library(dplyr)
library(knitr)
library(mia)
library(phyloseq)
# cargar el objeto tse y convertirlo a phyloseq
setwd("~/Dropbox/FCien/cursos/Metabarcoding_pedeciba")
tse <- readRDS("tse_dia4_base.rds")
tse
ps <- convertToPhyloseq(tse, assay.type = "counts")
ps
sample_data(ps)

#convertir el phyloseq a microeco
library(file2meco)
library(magrittr)
ps.eco <- phyloseq2meco(ps)
ps.eco$sample_table %<>% base::subset(Sampling.position == "Sticking" | 
                                        Sampling.position == "After singeing")

library(microeco)
library(tidyverse)
#preparar el formato de microeco y chequear los reads por muestra
dataset <- ps.eco 
dataset
dataset$tidy_dataset() 
print(dataset)
dataset$sample_sums()
dataset$sample_sums() %>% range #con esto veo la cantidad de secuencias por muestra
#y de paso verifico si estoy trabajando con la rarefaccionada o no

set.seed(123)
dataset$rarefy_samples(method = "rarefy", sample.size = 16525)
dataset$sample_sums() %>% range
dataset$cal_abund()

#Primer método de cálculo de abundancia diferencial####
diff.abund1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Sampling.position", 
                              taxa_level = "Genus", alpha = 0.05)

diff.abund1$plot_diff_bar(threshold = 4, group_order = c("Sticking", "After singeing"))
# we show 20 taxa with the highest LDA (log10)
diff.abund1$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("Sticking", "After singeing"))

g1 <- diff.abund1$plot_diff_bar(use_number = 1:20,group_order = c("Sticking", "After singeing"))
# plot the abundance using same taxa in g1
g2 <- diff.abund1$plot_diff_abund(select_taxa = diff.abund1$plot_diff_bar_taxa, plot_type = "barerrorbar"
                                  , add_sig = FALSE, errorbar_addpoint = FALSE, 
                                  errorbar_color_black = TRUE, group_order = c("Sticking", "After singeing"))
# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank())
p <- g1 %>% aplot::insert_right(g2)
p

#Segundo  método: ANOVA####
diff.abund2 <- trans_diff$new(dataset = dataset, method = "anova", group = "Sampling.position", 
                              taxa_level = "Genus", alpha = 0.05,
                              anova_post_test = "HSD.test")
diff.abund2$res_diff
write.csv(diff.abund2$res_diff, file = "diff_abund_anova.csv")
diff.abund2$plot_diff_abund(group_order = c("Sticking", "After singeing"))

#Tercer método: DESeq2####
diff.abund3 <- trans_diff$new(dataset = dataset, method = "DESeq2", group = "Sampling.position", 
                               alpha = 0.05)
diff.abund3$res_diff
write.csv(diff.abund3$res_diff, file = "diff_abund_DESeq2.csv")
