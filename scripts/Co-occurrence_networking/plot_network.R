library(ggplot2)
library(phyloseq)
library(ape)
library(microbiomeSeq)
library(igraph)
library(visNetwork)
library(tidygraph)
library(ggraph)

# Read in master data tables to create phyloseq object
ASV = otu_table(as.matrix(read.table("feature-table.tsv", check.names = F, sep="\t", row.names = 1, header = T, quote = "", stringsAsFactors = F)), taxa_are_rows=T)
TAX = tax_table(as.matrix(read.table("ASV_taxonomy.tsv", check.names = F, sep="\t", row.names = 1, header = T, quote = "", stringsAsFactors = F)))
MD = sample_data(read.table("hiseas_metadata.txt", check.names = F, sep="\t", header = T, row.names = 1, quote = "", stringsAsFactors = F))
physeq_master = merge_phyloseq(phyloseq(ASV,TAX,MD))
physeq_master_relative  = transform_sample_counts(physeq_master, function(x) x / sum(x) )

## Filter by a given parameter
physeq = subset_samples(physeq_master, feature == "composting toilet")   # indoor kitchen, bed base, living room, composting toilet
physeq = taxa_level(physeq, "ASV")
physeq = filter_taxa(physeq, function(x) var(x) > 0, TRUE)

# Calculate the co-occurrence network
co_occur = co_occurence_network(physeq, grouping_column = "feature", rhos = 0.7, method = "bicor", qval_threshold = 0.05, plotNetwork = F)

# igraph object not working: Error in toVisNetworkData(graph) : igraph must be a igraph object
# Export network to igraph object
graph = co_occur$net$graph
graph_data = toVisNetworkData(graph)

### Create aesthetics for nodes/edges 
# Extract flat tables from phyloseq object
taxonomy_table = as.data.frame(as(tax_table(physeq_master), "matrix"))
metadata = as.data.frame(as(sample_data(physeq_master), "matrix"))
count_table = as.data.frame(as(otu_table(physeq_master_relative), "matrix"))
# Calculate mean relabund for all ASVs across all samples
ASV_mean_abund = as.data.frame(rowMeans(count_table))
colnames(ASV_mean_abund) = "mean_relabund"
ASV_mean_abund$label = rownames(ASV_mean_abund)
# Get taxonomy
colnames(taxonomy_table)[8] = "label"
# Merge all aesthetics together
newnodes = graph_data$nodes
newnodes = merge(newnodes,taxonomy_table,by="label")
newnodes = merge(newnodes,ASV_mean_abund,by="label")

# After creates aesthetics, rename graph object
graph_data$nodes = newnodes

# Tidy up the graph nodes and edges into a table
routes_tidy = tbl_graph(nodes = graph_data$nodes, edges = graph_data$edges, directed = TRUE)

# Plot the network with the defined aesthetics
ASV_network = ggraph(routes_tidy, layout = "graphopt") + 
  geom_edge_link(aes(width = width)) + 
  geom_node_point(aes(
    size=size/2,
    fill=mean_relabund*100,
    alpha = 1,
    shape = 21,
    color="black"),
    stroke=1) +
  scale_edge_width(name="Pearson rho", range = c(0.25, 2)) +
  #scale_edge_alpha_identity() +
  #scale_edge_colour_identity() +
  scale_fill_gradient2(name="Rel Abund", low = "gray90",mid="gray80", high="black",midpoint = 0.01) +
  scale_color_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  scale_size_identity() +
  #geom_node_text(aes(label = class), size = 3, repel = F) +
  #labs(edge_width = "Letters") +
  theme_graph()

ASV_network



