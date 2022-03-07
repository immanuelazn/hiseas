#### EXPORTING FOR R ANALYSIS
#qiime tools export \
#--input-path frequency-filtered-surface-table-no-mito.qza \
#--output-path exported

#qiime tools export \
#--input-path taxonomy_290.qza \
#--output-path exported

#qiime tools export \
#--input-path tree_290.qza \
#--output-path exported

# Edit the column names exported/taxonomy.tsv
#sed \
#'1 s/Feature ID/#ASVID/;1 s/Taxon/taxonomy/;1 s/Confidence/confidence/' \
#exported/taxonomy.tsv \
#> exported/biom-taxonomy.tsv

# Combine taxonomy with BIOM data
#biom add-metadata \
#-i exported/feature-table.biom \
#-o exported/table-with-taxonomy.biom \
#--observation-metadata-fp exported/biom-taxonomy.tsv \
#--sc-separated taxonomy


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


# Calculate relative abundance
calculate_relative_abundance <- function(x) x / sum(x)

# Calculate geometric mean 
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Define the set of random numbers
set.seed(800)

# Define name of metadata file
metadata <- "hiseas_metadata.txt"
# test case
otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OTU = otu_table
TAX = tax_table(taxmat)
OTU


metadata  <- import_qiime_sample_data("hiseas_metadata.txt")
tree      <- read_tree_greengenes

# Read pathabun file
pathabun_file <- read.table("exported-resupply-test/pathabun-table.tsv",
                            header=TRUE, stringsAsFactors = FALSE, row.names=1)
# Remove x in front of column names due to column names starting with number
names(pathabun_file) <- sub("^X", "", names(pathabun_file))

# Create pretend taxa table 
pseudotax = matrix(sample(letters, 2772, replace = TRUE), nrow = nrow(pathabun_file), ncol = 7)
rownames(pseudotax) = rownames(pathabun_file)
colnames(pseudotax) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Export biom file and tree from QIIME2 and provide original metadata file
OTU = otu_table(pathabun_file, taxa_are_rows=TRUE)
TAX = tax_table(pseudotax)

# Combine all information into a single phyloseq object
physeq <- merge_phyloseq(OTU, TAX)


# Beta diversity PCoA plots
physeq_rar <- rarefy_even_depth(physeq, sample.size = 1694943)
bc_beta <- ordinate(physeq_rar, method = "PCoA", distance = "bray")


#Plot weighted unifrac PCoA
bc_plastic_wood <- plot_ordination(physeq_rar,
                bc_beta,
                color = "orig_env_material") +
  labs(title = "PCoA (weighted UniFrac)") +
  theme_bw() +
  stat_ellipse(type = "norm", size = 1) +
  scale_colour_manual(values = c("#80ccf2", "#fbb778"),
                      labels = c("Plastic", "Wood")) +
  guides(color = guide_legend("Surface Material")) +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), axis.title = element_text(size = 12), legend.text = element_text(size = 13), legend.title = element_text(size = 13))


ggsave(file="bc_plastic_wood.png", plot=bc_plastic_wood, 
         width=8, height=6)

