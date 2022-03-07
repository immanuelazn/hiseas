#R script to generate diversity, differential abundance and relative abundance plots

# Load CRAN packages
library(tidyverse)
library(vegan)
library(ape)

# Load Bioconductor packages
library(phyloseq)
library(DESeq2)

# Load additional ggplot packages
library(ggplot2)
library(ggthemes)
library(svglite)

# write function to calculate relative abundance
calculate_relative_abundance <- function(x) x / sum(x)

# write function to calculate geometric mean 
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# set seed for reproducibility
set.seed(800)

# import metadata file
metadata <- "hiseas_metadata.txt"

# Export biom file and tree from QIIME2 
# provide metadata file
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("hiseas_metadata.txt")
tree      <- read_tree_greengenes("tree.nwk")

# Convert tree from from multichotomous to dichotmous format
tree <- multi2di(tree)

# Combine information into a phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

# set readable column names
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",

#sanity check                                 "Genus", "Species")
head(tax_table(physeq))

# calculate sequencing depth, observe which samples have at least 14710 reads
sample_sums(physeq)
sample_sums(physeq) >= 14710

# Prune just to check
at_least_14710 <- prune_samples(sample_sums(physeq) >= 14710, physeq)
sample_sums(at_least_14710)

# Calculate relative abundance
physeq_RA <- transform_sample_counts(physeq, 
                                     calculate_relative_abundance)

# Total counts
total_counts <- taxa_sums(physeq)
relative_abundance <- calculate_relative_abundance(total_counts)

# Filter out low abundant features
abundant <- relative_abundance > 0.0005 
abundant_taxa <- prune_taxa(abundant, physeq)
abundant_taxa

#####differential abundance analysis#####

# subset by dust samples
dust <- subset_samples(physeq, env_material == "dust")

# find counts, remove low abundant features, and filter by metadata
dust_counts <- taxa_sums(dust)
relative_abundance_dust <- calculate_relative_abundance(dust_counts)
abundant_dust <- relative_abundance_dust > 0.0005
abundant_dust_taxa <- prune_taxa(abundant_dust, dust)

#set taxonomic analysis for the genus taxonomic rank
abundant_dust_genera <- tax_glom(abundant_dust_taxa, taxrank = "Genus")
abundant_dust_genera

# Differential abundance comparing the wood and plastic surfaces by defining as categorical and setting the reference groups
sample_data(abundant_dust_genera)$orig_env_material <-
  factor(sample_data(abundant_dust_genera)$orig_env_material,
         levels = c("wood", "plastic"))
 
# Create DESeq2 object, comparing the surface materials for differential abundance
deseq_dust <- phyloseq_to_deseq2(abundant_dust_genera, ~ orig_env_material)
geo_means <- apply(counts(deseq_dust), 1, calculate_gm_mean)
deseq_dust <- estimateSizeFactors(deseq_dust, geoMeans = geo_means)
deseq_dust <- DESeq(deseq_dust, fitType = "local")

dust_diff_abund <- results(deseq_dust)

# Set significance level, and filter for significant results only
alpha <- 0.05
significant_dust <- as.data.frame(dust_diff_abund)
significant_dust <- filter(significant_dust, padj < alpha)

# merge differentially abundant genera with taxonomic assignments, and sort by the highest log2FoldChange
genera_df <- as.data.frame(tax_table(abundant_dust_genera))
significant_dust <- merge(significant_dust, genera_df, by = "row.names")
significant_dust <- arrange(significant_dust, log2FoldChange)

#visualize data to check
dim(significant_dust)

#####visualizing differential abundance######

# convert "Genus" variable to a categorical "factor" variable for ggplot visualization
significant_dust <- mutate(significant_dust,
                          Genus = factor(Genus, levels = Genus))

# create bar plot for visualization of differential abundance between wood and plastic
wood_plastic_diff_abundance <- ggplot(significant_dust, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

# Save as svg for editing and enhancement
ggsave(file="wood_plastic_diff_abundance.svg", plot=wood_plastic_diff_abundance, width=10, height=10)

#####relative abundance analysis#####

# calculate relative abundance by normalizing within sample
dust_RA <- transform_sample_counts(dust, calculate_relative_abundance)

# remove low abundance features
dust_counts <- taxa_sums(dust)
relative_abundance_dust <- calculate_relative_abundance(dust_counts)
abundant_dust <- relative_abundance_dust > 0.0005
abundant_dust_RA_taxa <- prune_taxa(abundant_dust, dust_RA)

# agglomerate taxa information at genus level

abundant_dust_RA_genera <-  tax_glom(abundant_dust_RA_taxa, taxrank = "Genus")

# subset plastic sample genera found in differential abundance analysis for relative abundance analysis
Rothia <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Rothia")
Neisseria <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Neisseria")
Abiotrophia <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Abiotrophia")
Novosphingobium <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Novosphingobium")
Varibaculum <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Varibaculum")
Jeotgalicoccus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Jeotgalicoccus")
Eremococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Eremococcus")
Sphingopyxis <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Sphingopyxis")
Microbacterium <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Microbacterium")
Methylobacterium-Methylorubrum <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Methylobacterium-Methylorubrum")
Tsukamurella <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Tsukamurella")
Lactobacillus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Lactobacillus")


Fastidiosipila <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Fastidiosipila")
Granulicatella <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Granulicatella")
Porphyromonas <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Porphyromonas")
Prevotella <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Prevotella")
Megasphaera <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Megasphaera")
Peptoniphilus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Peptoniphilus")
Gemella <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Gemella")
Atopobium <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Atopobium")
Fenollaria <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Fenollaria")
Fusobacterium <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Fusobacterium")
Anaerococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Anaerococcus")
Ezakiella <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Ezakiella")
Gardnerella <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Gardnerella")
Cutibacterium <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Cutibacterium")
Finegoldia <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Finegoldia")
Staphylococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Staphylococcus")
Micrococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Micrococcus")
Enhydrobacter <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Enhydrobacter")
Acinetobacter <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Acinetobacter")
Corynebacterium <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Corynebacterium")
Ralstonia <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Ralstonia")

# subset wood sample genera found in differential abundance analysis for relative abundance analysis
Sphingobacterium <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Sphingobacterium")
Brevundimonas <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Brevundimonas")
Rhodococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Rhodococcus")
Enterococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Enterococcus")
Aerococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Aerococcus")
Massilia <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Massilia")
Duganella <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Duganella")
Chelatococcus <- subset_taxa(abundant_dust_RA_genera, Genus == "g__Chelatococcus")

#transform genera from wood samples to long format for ggplot
Sphingobacterium_long <- psmelt(Sphingobacterium)
Brevundimonas_long <- psmelt(Brevundimonas)
Rhodococcus_long <- psmelt(Rhodococcus)
Enterococcus_long <- psmelt(Enterococcus)
Aerococcus_long <- psmelt(Aerococcus)
Massilia_long <- psmelt(Massilia)
Duganella_long <- psmelt(Duganella)
Chelatococcus_long <- psmelt(Chelatococcus)

##visualize relative abundance of genera of interest from wood surfaces
Sphingobacterium_abundance_plot <- ggplot(Sphingobacterium_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c("#80ccf2", "#fbb778")) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13)) +
  labs(x     = "Surface material",
       y     = "Relative abundance")

ggsave(file="Sphingobacterium_abundance_plot.png", plot=Sphingobacterium_abundance_plot, 
       width=8, height=6)

Brevundimonas_abundance_plot <- ggplot(Brevundimonas_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Brevundimonas",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Brevundimonas_abundance_plot.svg", plot=Brevundimonas_abundance_plot, 
       width=8, height=6)

Rhodococcus_abundance_plot <- ggplot(Rhodococcus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Rhodococcus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Rhodococcus_abundance_plot.png", plot=Rhodococcus_abundance_plot, 
       width=8, height=6)

Enterococcus_abundance_plot <- ggplot(Enterococcus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Enterococcus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Enterococcus_abundance_plot.png", plot=Enterococcus_abundance_plot, 
       width=8, height=6)

Aerococcus_abundance_plot <- ggplot(Aerococcus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Aerococcus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Aerococcus_abundance_plot.png", plot=Aerococcus_abundance_plot, 
       width=8, height=6)

Massilia_abundance_plot <- ggplot(Massilia_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Massilia",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Massilia_abundance_plot.png", plot=Massilia_abundance_plot, 
       width=8, height=6)

Duganella_abundance_plot <- ggplot(Duganella_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Duganella",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Duganella_abundance_plot.png", plot=Duganella_abundance_plot, 
       width=8, height=6)

Chelatococcus_abundance_plot <- ggplot(Chelatococcus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Chelatococcus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Chelatococcus_abundance_plot.png", plot=Chelatococcus_abundance_plot, 
       width=8, height=6)


#transform genera from plastic samples to long format for ggplot
Rothia_long <- psmelt(Rothia)
Neisseria_long <- psmelt(Neisseria)
Abiotrophia_long <- psmelt(Abiotrophia)
Novosphingobium_long <- psmelt(Novosphingobium)
Varibaculum_long <- psmelt(Varibaculum)
Jeotgalicoccus_long <- psmelt(Jeotgalicoccus)
Eremococcus_long <- psmelt(Eremococcus)
Sphingopyxis_long <- psmelt(Sphingopyxis)
Microbacterium_long <- psmelt(Microbacterium)
Tsukamurella_long <- psmelt(Tsukamurella)
Lactobacillus_long <- psmelt(Lactobacillus)
Fastidiosipila_long <- psmelt(Fastidiosipila)
Granulicatella_long <- psmelt(Granulicatella)
Porphyromonas_long <- psmelt(Porphyromonas)
Prevotella_long <- psmelt(Prevotella)
Megasphaera_long <- psmelt(Megasphaera)
Peptoniphilus_long <- psmelt(Peptoniphilus)
Gemella_long <- psmelt(Gemella)
Atopobium_long <- psmelt(Atopobium)
Fenollaria_long <- psmelt(Fenollaria)
Fusobacterium_long <- psmelt(Fusobacterium)
Anaerococcus_long <- psmelt(Anaerococcus)
Ezakiella_long <- psmelt(Ezakiella)
Gardnerella_long <- psmelt(Gardnerella)
Cutibacterium_long <- psmelt(Cutibacterium)
Finegoldia_long <- psmelt(Finegoldia)
Staphylococcus_long <- psmelt(Staphylococcus)
Micrococcus_long <- psmelt(Micrococcus)
Enhydrobacter_long <- psmelt(Enhydrobacter)
Acinetobacter_long <- psmelt(Acinetobacter)
Corynebacterium_long <- psmelt(Corynebacterium)
Ralstonia_long <- psmelt(Ralstonia)

#visualize relative abundance of genera of interest from plastic surfaces
Fastidiosipila_abundance_plot <- ggplot(Fastidiosipila_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Fastidiosipila",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Fastidiosipila_abundance_plot.png", plot=Fastidiosipila_abundance_plot, 
       width=8, height=6)

Granulicatella_abundance_plot <- ggplot(Granulicatella_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Granulicatella",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Granulicatella_abundance_plot.png", plot=Granulicatella_abundance_plot, 
       width=8, height=6)

Porphyromonas_abundance_plot <- ggplot(Porphyromonas_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Porphyromonas",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Porphyromonas_abundance_plot.png", plot=Porphyromonas_abundance_plot, 
       width=8, height=6)

Prevotella_abundance_plot <- ggplot(Prevotella_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Prevotella",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Prevotella_abundance_plot.png", plot=Prevotella_abundance_plot, 
       width=8, height=6)

Megasphaera_abundance_plot <- ggplot(Megasphaera_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Megasphaera",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Megasphaera_abundance_plot.png", plot=Megasphaera_abundance_plot, 
       width=8, height=6)

Peptoniphilus_abundance_plot <- ggplot(Peptoniphilus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Peptoniphilus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Peptoniphilus_abundance_plot.png", plot=Peptoniphilus_abundance_plot, 
       width=8, height=6)

Gemella_abundance_plot <- ggplot(Gemella_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Gemella",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Gemella_abundance_plot.png", plot=Gemella_abundance_plot, 
       width=8, height=6)

Atopobium_abundance_plot <- ggplot(Atopobium_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Atopobium",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Atopobium_abundance_plot.png", plot=Atopobium_abundance_plot, 
       width=8, height=6)

Fenollaria_abundance_plot <- ggplot(Fenollaria_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Fenollaria",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Fenollaria_abundance_plot.png", plot=Fenollaria_abundance_plot, 
       width=8, height=6)

Fusobacterium_abundance_plot <- ggplot(Fusobacterium_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Fusobacterium",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Fusobacterium_abundance_plot.png", plot=Fusobacterium_abundance_plot, 
       width=8, height=6)

Anaerococcus_abundance_plot <- ggplot(Anaerococcus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Anaerococcus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Anaerococcus_abundance_plot.png", plot=Anaerococcus_abundance_plot, 
       width=8, height=6)

Ezakiella_abundance_plot <- ggplot(Ezakiella_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Ezakiella",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Ezakiella_abundance_plot.png", plot=Ezakiella_abundance_plot, 
       width=8, height=6)

Gardnerella_abundance_plot <- ggplot(Gardnerella_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Gardnerella",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Gardnerella_abundance_plot.png", plot=Gardnerella_abundance_plot, 
       width=8, height=6)

Cutibacterium_abundance_plot <- ggplot(Cutibacterium_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Cutibacterium",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Cutibacterium_abundance_plot.png", plot=Cutibacterium_abundance_plot, 
       width=8, height=6)
Finegoldia_abundance_plot <- ggplot(Finegoldia_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Finegoldia",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Finegoldia_abundance_plot.png", plot=Finegoldia_abundance_plot, 
       width=8, height=6)

Staphylococcus_abundance_plot <- ggplot(Staphylococcus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Staphylococcus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Staphylococcus_abundance_plot.svg", plot=Staphylococcus_abundance_plot, 
       width=8, height=6)

Micrococcus_abundance_plot <- ggplot(Micrococcus_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Micrococcus",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Micrococcus_abundance_plot.png", plot=Micrococcus_abundance_plot, 
       width=8, height=6)

Enhydrobacter_abundance_plot <- ggplot(Enhydrobacter_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Enhydrobacter",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Enhydrobacter_abundance_plot.png", plot=Enhydrobacter_abundance_plot, 
       width=8, height=6)

Acinetobacter_abundance_plot <- ggplot(Acinetobacter_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Acinetobacter",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Acinetobacter_abundance_plot.png", plot=Acinetobacter_abundance_plot, 
       width=8, height=6)

Corynebacterium_abundance_plot <- ggplot(Corynebacterium_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Corynebacterium",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Corynebacterium_abundance_plot.png", plot=Corynebacterium_abundance_plot, 
       width=8, height=6)

Ralstonia_abundance_plot <- ggplot(Ralstonia_long, aes(x = orig_env_material, y = Abundance, fill = orig_env_material)) +
  geom_boxplot() +
  labs(title = "Relative abundance of Ralstonia",
       x     = "Surface material",
       y     = "Relative abundance")
ggsave(file="Ralstonia_abundance_plot.png", plot=Ralstonia_abundance_plot, 
       width=8, height=6)

####end####
