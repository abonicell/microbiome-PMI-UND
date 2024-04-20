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

theme_set(theme_bw(16))

ps1 <- readRDS("ps1.rds")

# group taxa by similarity for both location
phy <- phyloseq::tax_glom(ps1, "Phylum")
(phyloseq::taxa_names(phy) <- phyloseq::tax_table(phy)[, "Phylum"])

# merge by week
time_merge <- merge_samples(phy, "week")

ASVnames10 = names(sort(taxa_sums(time_merge), TRUE)[1:10])
Time10  = prune_taxa(ASVnames10,  phy)
Time10 = prune_taxa(ASVnames10, time_merge)
print(Time10)

Time10@sam_data$week <- rownames(Time10@sam_data)

# transform in percentage
trans <-
  transform_sample_counts(Time10, function(x)
    (x / sum(x)) * 100)

# arrange PMI in increasing order
sample_data(trans)$week <-
  factor(
    sample_data(trans)$week,
    levels = c(
      "week_1",
      "week_2",
      "week_3",
      "week_4",
      "week_5",
      "week_6",
      "week_7",
      "week_8",
      "week_9",
      "week_10",
      "week_11",
      "week_12",
      "week_13",
      "week_14",
      "week_15",
      "week_16",
      "week_17",
      "week_18",
      "week_19",
      "week_20",
      "week_21",
      "week_22",
      "week_23"
    )
  )

bar_phy <- plot_bar(trans, x = "week", fill = "Phylum") +
  geom_bar(aes(fill = Phylum),
           stat = "identity",
           position = "stack") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(x = "PMI (weeks)", y = "Abundance (%)", title = "Total sample") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  )

bar_phy

# abundance
phyloseq::psmelt(trans) %>%
  group_by(OTU, week) %>% tally(Abundance) %>%
  write_csv("PMI_phylum.csv")

#-------------------------------------------------------------------------------
# group taxa by similarity
cls <- phyloseq::tax_glom(ps1, "Class")

# merge by week
time_merge <- merge_samples(cls, "week")

ASVnames10 = names(sort(taxa_sums(time_merge), TRUE)[1:10])
Time10  = prune_taxa(ASVnames10,  cls)
Time10 = prune_taxa(ASVnames10, time_merge)
print(Time10)

Time10@sam_data$week <- rownames(Time10@sam_data)

# transformed in percentage
trans <-
  transform_sample_counts(Time10, function(x)
    (x / sum(x)) * 100)

# arrange PMI in increasing order
sample_data(trans)$week <-
  factor(
    sample_data(trans)$week,
    levels = c(
      "week_1",
      "week_2",
      "week_3",
      "week_4",
      "week_5",
      "week_6",
      "week_7",
      "week_8",
      "week_9",
      "week_10",
      "week_11",
      "week_12",
      "week_13",
      "week_14",
      "week_15",
      "week_16",
      "week_17",
      "week_18",
      "week_19",
      "week_20",
      "week_21",
      "week_22",
      "week_23"
    )
  )

bar_class <- plot_bar(trans, x = "week", fill = "Class") +
  geom_bar(aes(fill = Class),
           stat = "identity",
           position = "stack") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(x = "PMI (weeks)", y = "Abundance (%)", title = "Total sample") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  )

bar_class

# abundance
phyloseq::psmelt(trans) %>%
  group_by(OTU, week) %>% tally(Abundance) %>%
  write_csv("time_class.csv")

#-------------------------------------------------------------------------------
# group taxa by similarity
# internal
ps1_int <- subset_samples(ps1, location == "Interior")

phy <- phyloseq::tax_glom(ps1_int, "Phylum")
(phyloseq::taxa_names(phy) <- phyloseq::tax_table(phy)[, "Phylum"])

# Time
time_merge <- merge_samples(phy, "week")

ASVnames10 = names(sort(taxa_sums(time_merge), TRUE)[1:10])
Time10  = prune_taxa(ASVnames10,  phy)
Time10 = prune_taxa(ASVnames10, time_merge)
print(Time10)

Time10@sam_data$week <- rownames(Time10@sam_data)

# transform in percentage values
trans <-
  transform_sample_counts(Time10, function(x)
    (x / sum(x)) * 100)

# arrange PMI in increasing order
sample_data(trans)$week <-
  factor(
    sample_data(trans)$week,
    levels = c(
      "week_1",
      "week_2",
      "week_3",
      "week_4",
      "week_5",
      "week_6",
      "week_7",
      "week_8",
      "week_9",
      "week_10",
      "week_11",
      "week_12",
      "week_13",
      "week_14",
      "week_15",
      "week_16",
      "week_17",
      "week_18",
      "week_19",
      "week_20",
      "week_21",
      "week_22",
      "week_23"
    )
  )

bar_phy_int <- plot_bar(trans, x = "week", fill = "Phylum") +
  geom_bar(aes(fill = Phylum),
           stat = "identity",
           position = "stack") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(x = "PMI (weeks)", y = "Abundance (%)", title = "Internal sample") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  )

bar_phy_int

# abundance
phyloseq::psmelt(trans) %>%
  group_by(OTU, week) %>% tally(Abundance) %>%
  write_csv("PMI_phylum_int.csv")

#-------------------------------------------------------------------------------
# group taxa by similarity
cls <- phyloseq::tax_glom(ps1_int, "Class")

# Time
Time_merge <- merge_samples(cls, "week")

ASVnames10 = names(sort(taxa_sums(Time_merge), TRUE)[1:10])
Time10  = prune_taxa(ASVnames10,  cls)
Time10 = prune_taxa(ASVnames10, Time_merge)
print(Time10)

Time10@sam_data$week <- rownames(Time10@sam_data)

trans <-
  transform_sample_counts(Time10, function(x)
    (x / sum(x)) * 100)

# arrange PMI in increasing order
sample_data(trans)$week <-
  factor(
    sample_data(trans)$week,
    levels = c(
      "week_1",
      "week_2",
      "week_3",
      "week_4",
      "week_5",
      "week_6",
      "week_7",
      "week_8",
      "week_9",
      "week_10",
      "week_11",
      "week_12",
      "week_13",
      "week_14",
      "week_15",
      "week_16",
      "week_17",
      "week_18",
      "week_19",
      "week_20",
      "week_21",
      "week_22",
      "week_23"
    )
  )

bar_class_int <- plot_bar(trans, x = "week", fill = "Class") +
  geom_bar(aes(fill = Class),
           stat = "identity",
           position = "stack") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(x = "PMI (hrs)", y = "Abundance (%)", title = "Internal sample") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  )

bar_class_int

# abundance
phyloseq::psmelt(trans) %>%
  group_by(OTU, week) %>% tally(Abundance) %>%
  write_csv("Time_class_int.csv")

#-------------------------------------------------------------------------------
# group taxa by similarity
# internal
ps1_ext <- subset_samples(ps1, location == "Exterior")

phy <- phyloseq::tax_glom(ps1_ext, "Phylum")
(phyloseq::taxa_names(phy) <- phyloseq::tax_table(phy)[, "Phylum"])

# Time
Time_merge <- merge_samples(phy, "week")

ASVnames10 = names(sort(taxa_sums(Time_merge), TRUE)[1:10])
Time10  = prune_taxa(ASVnames10,  phy)
Time10 = prune_taxa(ASVnames10, Time_merge)
print(Time10)

Time10@sam_data$week <- rownames(Time10@sam_data)

trans <-
  transform_sample_counts(Time10, function(x)
    (x / sum(x)) * 100)

# arrange PMI in increasing order
sample_data(trans)$week <-
  factor(
    sample_data(trans)$week,
    levels = c(
      "week_1",
      "week_2",
      "week_3",
      "week_4",
      "week_5",
      "week_6",
      "week_7",
      "week_8",
      "week_9",
      "week_10",
      "week_11",
      "week_12",
      "week_13",
      "week_14",
      "week_15",
      "week_16",
      "week_17",
      "week_18",
      "week_19",
      "week_20",
      "week_21",
      "week_22",
      "week_23"
    )
  )

bar_phy_ext <- plot_bar(trans, x = "week", fill = "Phylum") +
  geom_bar(aes(fill = Phylum),
           stat = "identity",
           position = "stack") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(x = "PMI (hrs)", y = "Abundance (%)", title = "External sample") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  )

bar_phy_ext

# abundance
phyloseq::psmelt(trans) %>%
  group_by(OTU, week) %>% tally(Abundance) %>%
  write_csv("PMI_phylum_ext.csv")

###############################################################################
# group taxa by similarity
cls <- phyloseq::tax_glom(ps1_ext, "Class")

# Time
Time_merge <- merge_samples(cls, "week")

ASVnames10 = names(sort(taxa_sums(Time_merge), TRUE)[1:10])
Time10  = prune_taxa(ASVnames10,  cls)
Time10 = prune_taxa(ASVnames10, Time_merge)
print(Time10)

Time10@sam_data$week <- rownames(Time10@sam_data)

trans <-
  transform_sample_counts(Time10, function(x)
    (x / sum(x)) * 100)

# arrange PMI in increasing order
sample_data(trans)$week <-
  factor(
    sample_data(trans)$week,
    levels = c(
      "week_1",
      "week_2",
      "week_3",
      "week_4",
      "week_5",
      "week_6",
      "week_7",
      "week_8",
      "week_9",
      "week_10",
      "week_11",
      "week_12",
      "week_13",
      "week_14",
      "week_15",
      "week_16",
      "week_17",
      "week_18",
      "week_19",
      "week_20",
      "week_21",
      "week_22",
      "week_23"
    )
  )

bar_class_ext <- plot_bar(trans, x = "week", fill = "Class") +
  geom_bar(aes(fill = Class),
           stat = "identity",
           position = "stack") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(x = "PMI (hrs)", y = "Abundance (%)", title = "External sample") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 0.9
    )
  )

bar_class_ext

# abundance
phyloseq::psmelt(trans) %>%
  group_by(OTU, week) %>% tally(Abundance) %>%
  write_csv("Time_class_ext.csv")

# A, C, and D reports phyla and the remaining column of plot class
ggarrange(
  bar_phy,
  bar_class,
  bar_phy_int,
  bar_class_int,
  bar_phy_ext,
  bar_class_ext,
  labels="AUTO",
  nrow = 3,
  ncol = 2
)

ggsave("barplots.pdf", width = 16, height = 16)

