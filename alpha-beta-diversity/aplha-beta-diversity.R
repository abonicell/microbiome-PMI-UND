# load packages
suppressPackageStartupMessages({
library("qiime2R") 
library("phyloseq")
library("tidyverse")
library("ggrepel")
library("RColorBrewer")
library("ggpubr")
library("vegan")
library("ranacapa")
library("caret")
library("knitr")
library("doParallel")
library("mlbench")
library("randomForestExplainer")
library("ranger")
library("microbiome")
library("microViz")
library("tidymodels")
library("indicspecies")
library("pheatmap")
library("ggplotify")
library("patchwork")
  
})

theme_set(theme_bw(14))

ps1 <- readRDS("ps1.rds")

# group taxa by similarity
phy <- phyloseq::tax_glom(ps1, "Phylum")
( phyloseq::taxa_names(phy) <- phyloseq::tax_table(phy)[, "Phylum"] )

# boxplot of taxa (phylum level)
phyloseq::psmelt(phy) %>%
  ggplot(data = ., aes(x = location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "Sampling Location", y = "Abundance\n") +
  facet_wrap( ~ OTU, scales = "free_y") +
  labs(title = NULL, subtitle = NULL) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) + theme_bw() +
  theme(legend.position = 'bottom')+
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1))

ggsave('phylum_box_Location.pdf', width = 11, height = 9)

# boxplot of taxa (phylum level)
phyloseq::psmelt(phy) %>%
  ggplot(data = ., aes(x = week, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "PMI (Times)", y = "Abundance\n") +
  facet_wrap( ~ OTU, scales = "free_y") +
  labs(title = NULL, subtitle = NULL) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) + theme_bw() +
  theme(legend.position = 'bottom')+
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(size = rel(0.8),angle=45,hjust=1,vjust=1.1))

ggsave('phylum_box_Time.pdf', width = 20, height = 10)

phyloseq::psmelt(phy) %>%
  ggplot(data = ., aes(x = Pig, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "Sampling Location", y = "Abundance\n") +
  facet_wrap( ~ OTU, scales = "free_y") +
  labs(title = NULL, subtitle = NULL) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) + 
  theme_bw() +
  theme(legend.position = 'bottom')+
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,vjust=0.9))

ggsave('phylum_box_pig.pdf', width = 10, height = 9)

#-------------------------------------------------------------------------------
# calculate alpha diversity - variation within sample
# both Location combined according to pig

my_comparisons <- list( c("P1","P2"),c("P1","P2"),c("P2","P3"))

rich_pig <- plot_richness(
  ps1,
  x = "pig",
  color = "pig",
  scales = "free_y",
  measures = c("Observed", "Shannon")
) +
  geom_boxplot() +
  ggpubr::stat_compare_means(label = "p.signif",method="t.test",
                             comparisons = my_comparisons,label.x = 1)+
  theme_bw() +
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)

rich_pig

my_comparisons <- list( c("Internal","External"))

# both Location combined according to Location
rich_Location <- plot_richness(
  ps1,
  x = "location",
  color = "location",
  scales = "free_y",
  measures = c("Observed", "Shannon")
) +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = my_comparisons,label = "p.signif", method="t.test", label.x = 1) +
  theme_bw()+
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)

rich_Location

# both Location combined according to Location
rich_Time <- plot_richness(
  ps1,
  x = "week",
  color = "week",
  scales = "free_y",
  measures = c("Observed", "Shannon")
) +
  geom_boxplot() +
  theme_bw() + theme(axis.text.x = element_text(
    size = rel(0.8),
    angle = 45,
    hjust = 1
  )) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)

rich_Time

#------------------------------------------------------------------------------
# unifrac distance as test - beta div = between group
unwunifrac_dist = phyloseq::distance(ps1, method = "unifrac", weighted = F)
ordination = ordinate(ps1, method = "PCoA", distance = unwunifrac_dist)

uni_pig <- plot_ordination(ps1, ordination, color = "pig")   +
  geom_point(aes(color = `pig`), alpha = 0.5, size = 4)   +
  theme_bw()+
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)

uni_pig

ggarrange(rich_pig, uni_pig,
          labels = 'AUTO')

ggsave("pig_eval.pdf", width = 14.10, height = 5.33)

uni_Location <- plot_ordination(ps1, ordination, color = "location")   +
  geom_point(aes(color = `location`), alpha = 0.5, size = 4)   +
  theme_bw()+
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)
  
plot_ordination(ps1, ordination, color = "week")   +
  geom_point(aes(color = `week`), alpha = 0.5, size = 4)   +
  theme_bw() +
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE) -> uni_Time

uni_Time

ggarrange(rich_Location, uni_Location,
          labels = 'AUTO')

ggsave("Location_eval.pdf", width = 14.10, height = 5.33)

ggarrange(rich_Time, uni_Time,
          labels = 'AUTO')

ggsave("Time_eval.pdf", width = 20.10, height = 5.33)
