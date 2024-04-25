# load packages
suppressPackageStartupMessages({
library("qiime2R") 
library("phyloseq")
library("tidyverse")
library("ggpubr")
library("vegan")
library("ranacapa")
library("microViz")
})

theme_set(theme_bw(16))

# read other files and merge into phyloseq object
taxonomy<-qiime2R::read_qza("taxonomy-silva.qza")
head(taxonomy$data);
qiime2R::parse_taxonomy(taxonomy$data)

( metadata <- readr::read_tsv("metadata.txt") )

( physeq<-qiime2R::qza_to_phyloseq(
  features="table-OTU97.qza",
  tree = "rooted-tree.qza",
  taxonomy="taxonomy-silva.qza",
  metadata = "metadata.txt"
) )# 4285 taxa

# remove contaminants
( physeq <- physeq %>%
    subset_taxa(
      Family!= "Mitochondria" | is.na(Family) #&
      # Order!="Chloroplast" | is.na(Class)
    ) %>%
    subset_taxa(
      # Family!= "Mitochondria" | is.na(Family) &
      Order!="Chloroplast" | is.na(Class)
    ) )# 4173 taxa

# prune sample to remove all with the sum of 0
ps_prune <-  prune_samples(sample_sums(physeq)!=0, physeq)
ps_prune # removes 8 samples with no values

# check for missing data
any(is.na(ps_prune@otu_table))

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(ps_prune))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps_prune, standf)

# check for missing data
any(is.na(ps@otu_table))

# arrange PMI in increasing order
sample_data(ps)$week <-
  factor(
    sample_data(ps)$week,
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

# replace unique identifier with consecutive ASV name
ps1 <-microViz::tax_names2rank(ps, colname = "unique")
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))

# save object for further processing in each subforlder for further porcessing
saveRDS(ps1, file = "/Users/andreabonicelli/Documents/Microbiome/lavina_pigs/microbiome_UND/bar-plots/ps1.rds")
saveRDS(ps1, file = "/Users/andreabonicelli/Documents/Microbiome/lavina_pigs/microbiome_UND/alpha-beta-diversity/ps1.rds")
saveRDS(ps1, file = "/Users/andreabonicelli/Documents/Microbiome/lavina_pigs/microbiome_UND/estimation-models/ps1.rds")

