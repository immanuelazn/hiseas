# Load CRAN packages
library(tidyverse)
library(vegan)
library(ape)

# Load Bioconductor packages

library(DESeq2)
library(biomformat)
library(phyloseq)


# Load additional ggplot packages
library(ggplot2)
library(ggthemes)
library(svglite)


# Define the set of random numbers
set.seed(820)

###################################################################################
# Helper Function given by Gavin Douglas for Level 3 KO categorization
categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.
  
  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))
  
  colnames(out_pathway) <- c("pathway", colnames(in_ko))
  
  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}
#######################################################################################
# Count how many ASVs filtered

# Read raw table of NSTIs
nsti_asv <- read.table("marker_predicted_and_nsti.tsv", header = TRUE)
# Determine how many sequences there are
nrow(nsti_asv)
# 14021 total asvs
nrow(nsti_asv[nsti_asv$metadata_NSTI  > 2,])
# 228 asvs that were filtered out

### Count how many samples have NSTI Values below 

# Count how many samples have below 0.15 NSTI
nsti_sample <- read.table("weighted_nsti.tsv", header = TRUE)
nrow(nsti_sample)
#104 total samples
nrow(nsti_sample[nsti_sample$weighted_NSTI < 0.15,])
#104 samples with NSTI under 104.
#######################################################################################
# Convert kegg with descrips from tsv to tsv for STAMP.  Remove KO and only keep descrip

# Read ko tsv that had commas and spaces replaced
ko_descript <- read.table("pred_metagenome_unstrat_descrip.tsv",
                             header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)
# Remove x in front of column names due to column names starting with number
names(ko_descript) <- sub("^X", "", names(ko_descript))

# Turn table into matrix
otumat_ko = as(ko_descript, "matrix")
# Remove error duplicate KO
otumat_ko<- otumat_ko[-2298,]

# Convert descript column into row names
otumat_ko_row <- otumat_ko[,1]
rownames(otumat_ko) <- otumat_ko_row

#Remove descript  column
otumat_ko <- otumat_ko[,-1]
# Convert to numeric matrix
otumat_ko_numeric <- mapply(otumat_ko, FUN=as.numeric)
otumat_ko_new <- matrix(data=otumat_ko_numeric, ncol =104, nrow = 6481)
rownames(otumat_ko_new) <- rownames(otumat_ko)
colnames(otumat_ko_new) <- colnames(otumat_ko)


# export into tsv
write.table(otu_table(otumat_ko_new, taxa_are_rows=TRUE), "kegg_descrip.tsv",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Go to tsv and add "descriptions" on line 1"

#####################################################################################
# Convert pathabun with descrips from csv to tsv for STAMP.  Remove pwy and only keep descrip

# Read pathway csv that had commas and spaces replaced
abund_descript <- read.csv("pathway-abund-descrip-comma.csv", header = TRUE)
# Remove x in front of column names due to column names starting with number
names(abund_descript) <- sub("^X", "", names(abund_descript))

# Turn table into matrix
otumat_abund = as(abund_descript, "matrix")

# Remove pwy column
otumat_abund <- otumat_abund[, -1]

# Convert processes column into row names
otumat_abund_row <- otumat_abund[,1]
rownames(otumat_abund) <- otumat_abund_row

#Remove process column
otumat_abund <- otumat_abund [,-1]
# Convert to numeric matrix
otumat_abund_numeric <- mapply(otumat_abund, FUN=as.numeric)
otumat_abund_new <- matrix(data=otumat_abund_numeric, ncol =104, nrow = 392)
rownames(otumat_abund_new) <- rownames(otumat_abund)
colnames(otumat_abund_new) <- colnames(otumat_abund)


# export into tsv
write.table(otu_table(otumat_abund_new, taxa_are_rows=TRUE), "pathabun_descrip.tsv",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Go to tsv and add "descriptions" on line 1"


######################################################
#################################################################################
# KO Level 3 BRITE Mapping


# Import KO PICRUSt outputs and KEGG mapping
ko_data <- read.table("pred_metagenome_unstrat.tsv", header=TRUE, sep="\t", row.names=1)

kegg_brite_map <- read.table("picrust1_KO_BRITE_map.tsv",
                             header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)

#Run function to categorize all KOs by level 3 in BRITE hierarchy.
ko_data_l3 <- categorize_by_function_l3(ko_data, kegg_brite_map)

# Remove x in front of column names due to inconsistent formatting
names(ko_data_l3) <- sub("^X", "", names(ko_data_l3))

# Turn table into matrix
otumat_kegg = as(ko_data_l3, "matrix")

# Import data for phyloseq
OTU = otu_table(otumat_kegg, taxa_are_rows=TRUE)
metadata  <- import_qiime_sample_data("hiseas_metadata.txt")
phylo <- merge_phyloseq(OTU, metadata)

# Potential DESEq analysis here once I get it working

# Export ko otu-table for STAMP
otu_export<-as(otu_table(phylo),"matrix")
otu_biom<-make_biom(data=otu_export)
write.table(OTU, "l3-ko.tsv",sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
### Go into the tsv and add a "pathway" to col names on line 1


