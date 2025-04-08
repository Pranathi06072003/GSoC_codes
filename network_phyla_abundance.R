
##merge_data1 taken from SPIEC EASI Function end code 
ps_filtered <- merge_data1
otu.table <- as(otu_table(ps_filtered), "matrix")
env.all <- as(sample_data(ps_filtered), "data.frame")
otu.table.all <- t(otu.table)

#row.names(env.all) == row.names(otu.table.all)

if (nrow(otu.table) != nsamples(merge_data)) {
  # If samples are not rows, transpose the OTU table
  otu.table <- t(otu.table)
}
otu.names <- colnames(otu.table)
num_cores <- parallel::detectCores() - 1  # Use one less than the total number of cores to avoid system overload

run_spiec_easi <- function(otu_table, num_cores) {
  future.apply::future_lapply(1, function(x) {
    spiec.easi(
      otu_table,
      method = 'mb',
      lambda.min.ratio = 1e-2,
      nlambda = 20,
      icov.select.params = list(rep.num = 50, ncores = 1)  # ncores = 1 because future will handle parallelization
    )
  }, future.seed = TRUE)[[1]]  # Set future.seed = TRUE as an argument to future_lapply
}

# Set up the parallel backend
spieceasi.net <-  run_spiec_easi(otu.table, num_cores)
# The result is an object of class "spiec.easi"
print(spieceasi.net)

# Reset the plan to default
plan(sequential)

adjacency_matrix <- symBeta(getOptBeta(spieceasi.net))
colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- colnames(otu.table)

vsize <- log2(apply(otu.table, 2, mean)) # add log abundance as properties of vertex/nodes.

marine.ig <- graph.adjacency(adjacency_matrix, mode='undirected', add.rownames = TRUE, weighted = TRUE)
otu.names <- colnames(otu.table)
V(marine.ig)$name <- otu.names

# Extract taxonomy information

taxa <- tax_table(ps_filtered)
phylum <- as.factor(taxa[, "Phylum"])
phylum <- gsub("^p__", "", phylum)

# Convert phylum to a factor
phylum <- factor(phylum)

# Consolidate phylum counts
phylum_counts <- table(phylum)

# Consolidate repeated phylum names
unique_phylum <- unique(phylum)
phylum_levels <- levels(unique_phylum)
consolidated_counts <- sapply(phylum_levels, function(level) sum(phylum_counts[names(phylum_counts) == level]))

# Calculate the proportions of each phylum in the network
phylum_proportions <- prop.table(consolidated_counts) * 100

# Print the consolidated counts and proportions
print(consolidated_counts)
print(phylum_proportions)


unique_phyla <- levels(phylum)
color_palette <- brewer.pal(n = min(length(unique_phyla), 12), name = "Paired")
if(length(unique_phyla) > length(color_palette)) {
  color_palette <- colorRampPalette(color_palette)(length(unique_phyla))
}
taxa_colors <- c(color_palette, "grey")
names(taxa_colors) <- unique_phyla

par(mar = c(1, 1, 1, 1) + 0.1)  # Increase the right margin for the legend
#layout <- layout_with_fr(marine.ig)
plot(marine.ig,
     vertex.size=5, 
     vertex.label=NA, 
     edge.color="grey",
     vertex.color=taxa_colors[phylum],
     main="Network Colored by Taxonomy")
legend("topright", legend=names(taxa_colors), 
       col=taxa_colors, pch=19, pt.cex=2, cex=0.8, bty="n", inset=c(-0.1, 0))
#dev.off()

# Add legend for taxonomy colors
legend("topright", legend=levels(phylum), 
       col=taxa_colors, pch=19, pt.cex=2, cex=0.8, bty="n")


# Define colors for the plot
taxa_colors <- brewer.pal(n = 30, name = "Paired")
par(mar = c(10, 5, 4, 0) + 0.1) 

#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/NetworkProperties/Phylum_Barplot.png", width = 800, height = 600)

# Visualize the proportions using a bar plot
barplot(phylum_proportions, col=taxa_colors, las=2,
        main="Proportion of Different Phyla in the Network",
        ylab="Percentage")
#dev.off()
# Optionally, visualize the proportions using a pie chart

par(mar = c(1, 0, 2, 2) + 0.1) 
pie(phylum_proportions, labels=names(phylum_proportions),
    col=taxa_colors, main="Proportion of Different Phyla in the Network")

# Plot the network colored by taxonomy
adjacency_matrix <- as.matrix(adjacency_matrix)
adjacency_matrix.m <- melt(adjacency_matrix)
node_taxonomy <- unique(c(as.character(unique(adjacency_matrix.m$Var1)),as.character(unique(adjacency_matrix.m$Var2))))
node_ids <- c(1:length(node_taxonomy))
node_labels <- c(1:length(node_taxonomy))

phylum_vector <- as.character(phylum[match(node_ids, V(marine.ig)$name)])
nodes_phylum <- phylum_vector

node_data <- data.frame(
  ID = node_ids,
  Label = node_labels,
  Taxonomy = node_taxonomy,
  Phylum = phylum_vector
)
write.csv(node_data, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Nodes_538Nodes.csv")
