# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)  # For performing PCA
library(dplyr)
#library(NetCoMi)

remove_prefixes <- function(physeq_1) {
  tax_table1 = tax_table(physeq_1)
  taxa_names <- gsub("^k_", "", tax_table1) # Remove 'k_' prefix
  taxa_names <- gsub("^p_", "", taxa_names)
  taxa_names <- gsub("^c_", "", taxa_names)
  taxa_names <- gsub("^o_", "", taxa_names)
  taxa_names <- gsub("^f_", "", taxa_names)
  taxa_names <- gsub("^g_", "", taxa_names)
  taxa_names <- gsub("^s_", "", taxa_names) # Remove 's_' prefix
  taxa_names <- gsub("^_", "", taxa_names) # Remove 's_' prefix
  tax_names2 <- tax_table(taxa_names)
  otu_table1 <- otu_table(physeq_1)
  metadata <- sam_data(physeq_1)
  physeq2 <- phyloseq(otu_table1, tax_names2, metadata)
  
  
  return(physeq2)
}

combined_data1 <- remove_prefixes(combined_data)
combined_genus_agglomeration1 <- tax_glom(combined_data1, taxrank = 'Genus')
#combined_genus_agglomeration <- tax_glom(combined_data, taxrank = 'Genus')


complete_data1 <- prune_taxa(taxa_sums(combined_genus_agglomeration1) > 0, combined_genus_agglomeration1)

#complete_data1 <- core(complete_data1, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data1 <- prune_samples(non_zero_samples, complete_data1)



combined_N <- transform_sample_counts(complete_data1, function(x) x / sum(x))

#merged_NMDS <- ordinate(combined_data1, method = "NMDS", distance = "bray")
merged_PCoA <- ordinate(complete_data1, method = "PCoA", distance = "bray")
merged_PCoA_N <- ordinate(combined_N, method = "PCoA", distance = "bray")

plot_ordination(combined_N, merged_NMDS, type = "samples", color = "sample_origin")
plot_ordination(complete_data1, merged_PCoA, axes= c(1,2), type= "samples", color="sample_origin", title = "Principal Component Analysis")
plot_ordination(combined_N, merged_PCoA_N, axes= c(1,2), type= "samples", color="sample_origin", title = "Principal Component Analysis After Normalisation")
plot_ordination(combined_N, merged_PCoA_N, axes= c(1,3), type= "samples", color="sample_origin", title = "Principal Component Analysis After Normalisation")

eigen_values <- merged_PCoA_N$values$Eigenvalues
eigenvalues_perc <- eigen_values/sum(eigen_values)

pca_df <- data.frame(
  Principal_Component <- factor(pca_df$Principal_Component, levels = paste0("PC", 1:10)), 
  Proportion_Variance_Explained <- eigenvalues_perc[1:10]*100  # in percentage
)

#pca_df$Principal_Component <- factor(pca_df$Principal_Component, levels = 1:10)


ggplot(pca_df, aes(x = Principal_Component)) +
  geom_bar(aes(y = Proportion_Variance_Explained), stat = "identity", fill = "violet", alpha = 0.7) +
  labs(
    title = "Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Percentage of Variance Explained"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))  # Center the title
