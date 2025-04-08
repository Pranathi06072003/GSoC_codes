library(SpiecEasi)
library(SpiecEasi)
library(devtools)
library(igraph)
library(vegan)
library(Matrix)
library(reshape2)
library(plyr)
library(future)
library(future.apply)
library(dplyr)
library(gridExtra)
library(grid)
library(igraph)
library(intergraph)
library(GGally)
library(network)
library(intergraph)
library(RColorBrewer)
library(dplyr)
library(microbiome)
library(NetCoMi)

node_gen <- function(merge_data1, filepath_node, filepath_edge, metrics_path, degree_path){
  
  merge_data1.f <- microbiomeutilities::format_to_besthit(merge_data1)
  otu.table <- as(otu_table(merge_data1.f), "matrix")
  env.all <- as(sample_data(merge_data1.f), "data.frame")
  otu.table.all <- t(otu.table)
  row.names(env.all) == row.names(otu.table.all)
  
  if (nrow(otu.table) != nsamples(merge_data1.f)) {
    # If samples are not rows, transpose the OTU table
    otu.table <- t(otu.table)
    
  }
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
  
  ##Generating the adjacency matrix
  adjacency_matrix <- symBeta(getOptBeta(spieceasi.net))
  colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- colnames(otu.table)
  adjacency_matrix_return <- as.matrix(adjacency_matrix)
  
  vsize <- log2(apply(otu.table, 2, mean)) # add log abundance as properties of vertex/nodes.
  
  
  marine.ig <- graph.adjacency(adjacency_matrix, mode='undirected', add.rownames = TRUE, weighted = TRUE)
  otu.names <- colnames(otu.table)
  V(marine.ig)$name <- otu.names
  ##generates all the attributes and wieghts
  
  marine.net <- asNetwork(marine.ig)
  network::set.edge.attribute(marine.net, "color", ifelse(marine.net %e% "weight" > 0, "steelblue", "orange"))
  
  
  phyla <- map_levels(colnames(otu.table), from = "best_hit", to = "Phylum", tax_table(merge_data1.f))
  marine.net %v% "Phylum" <- phyla
  marine.net %v% "nodesize" <- vsize
  
  ##converting the weights to dstances 
  net.dist <- marine.ig
  max(abs(E(net.dist)$weight))
  weights.dist <- 1 - abs(E(net.dist)$weight)
  E(net.dist)$weight <- weights.dist
  
  net.abs <- marine.ig
  E(net.abs)$weight <- abs(E(net.abs)$weight)
  
  ##modules
  wt <- cluster_louvain(marine.ig, weights = E(net.dist)$weight)
  temp <- V(marine.ig)$name
  temp <- as.data.frame(temp)
  temp$louvain <- membership(wt)
  V(marine.ig)$louvain <- temp$louvain
  
  length(unique(temp$louvain))
  summary_modules <- data.frame(table(temp$louvain))
  colnames(summary_modules) <- c("louvain", "n")
  summary_modules
  modules <- as.numeric(summary_modules$louvain[which(summary_modules$n >3)])
  
  x <- max(modules)+1
  for (i in c(1:length(temp$temp))) {
    if(!(temp$louvain[i] %in% modules)){
      temp$louvain[i] <- paste(x)
    }}
  modules <- temp
  modules$louvain <- as.numeric(modules$louvain)
  modules <- modules[order(modules$louvain),]
  module.lookup <- data.frame("louvain"=unique(modules$louvain),"new_louvain" = c(1:length(unique(modules$louvain))))
  new <- merge(modules,module.lookup)
  modules <- new
  modules <- modules[,2:3]
  summary_modules <- data.frame(table(modules$new_louvain))
  summary_modules
  max(modules$new_louvain)
  
  ##centrality measures
  # alpha centrality
  net.alpha <- alpha.centrality(marine.ig)
  # degree distribution
  net.strength <- strength(net.abs)
  # betweenness centrality
  bet <- betweenness(net.dist,v = V(net.dist))
  
  # make a summary of centrality metrics
  
  summary_cent <- as.data.frame(net.alpha)
  colnames(summary_cent) <- ("Alpha_centrality")
  rownames(summary_cent) <- colnames(otu.table)
  summary_cent$Weighted_vertex_degree <- net.strength
  summary_cent$Betweenness_centrality <- bet
  metrics <- summary_cent
  
  write.csv(metrics, metrics_path)
  
  adjacency_matrix <- as.matrix(adjacency_matrix)
  adjacency_matrix.m <- melt(adjacency_matrix)
  
  stl.mb <- degree.distribution(marine.ig)
  png(degree_path, width = 800, height = 600, res = 120)
  
  plot(0:(length(stl.mb)-1), stl.mb, ylim=c(0,.35), type='b', ylab="Frequency", xlab="Degree", main="Degree Distributions")
  
  colnames(adjacency_matrix.m) <- c("Var1", "Var2", "value")
  
  
  dev.off()
  
  node.names <- unique(c(as.character(unique(adjacency_matrix.m$Var1)),as.character(unique(adjacency_matrix.m$Var2))))
  node.names <- as.data.frame(node.names)
  node.names$node_number <- c(1:length(node.names$node.names))
  node.names$node_number2 <- c(1:length(node.names$node.names))
  colnames(node.names) <- c("Taxonomy", "Label", "Id")
  row.names(node.names) <- node.names$Taxonomy
  
  node.names <- node.names[, c("Taxonomy", "Label", "Id")]
  
  
  row.names(modules) <- modules$temp
  modules <- modules[order(modules$temp),]
  row.names(node.names) ==row.names(metrics)
  row.names(node.names) ==row.names(modules)
  node.names.final <- cbind(node.names, metrics,modules)
  
  node.names.final <- cbind(node.names, metrics,modules)
  write.table(node.names.final, filepath_node, sep = ",", row.names= FALSE)
  
  temp <- merge(x = adjacency_matrix.m, y = node.names, by.x = "Var1", by.y = "Taxonomy")
  # create the edge list
  colnames(temp) <- c("Var1","Var2","value","remove","source_number")
  temp <- temp[,-4]
  edge.list <- merge(x = temp, y = node.names, by.x = "Var2", by.y = "Taxonomy")
  colnames(edge.list) <- c("Var1","Var2","value","source.number","target.number")
  edge.list <- edge.list[,c(3,4,6)]
  colnames(edge.list) <- c("value","Var1","Var2")
  edge.list$Type <- "Undirected"
  negative <- ifelse(edge.list$value<0, "negative", "positive")
  edge.list$Negative <- negative
  edge.list$value <- abs(edge.list$value)
  edge.list <- edge.list[which(abs(edge.list$value)>0),]
  write.table(edge.list, filepath_edge, sep = ",", row.names = FALSE)
  
  return(adjacency_matrix_return)
}

##10% prevalence 


#complete_data1 <- prune_taxa(taxa_sums(combined_data) > 500, combined_data)
#complete_data1 <- tax_glom(complete_data1, taxrank = "Genus")

#combined_genus_agglomeration1 <- tax_glom(combined_data, taxrank = 'Genus')
complete_data1 <- core(combined_genus_agglomeration1, detection=0, prevalence = 0.05, include.lowest = FALSE)
non_zero_samples <- sample_sums(complete_data1) > 0
complete_data2 <- prune_samples(non_zero_samples, complete_data1)
#complete_data2 <- renameTaxa(complete_data2, 
                             #pat = "<name>", 
                             #substPat = "<subst_name>(<subst_R>)",
                             #numDupli = "Genus")

adj_0.05_Prev <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/degree_0.01_Prev.png")
net_0.05_Prev <- netConstruct(data=adj_0.01_Prev,
                             normMethod = "none", zeroMethod = "none",
                             sparsMethod = "none", dataType = "condDependence",
                             verbose = 3)
combined_analysis_0.05 <- netAnalyze(net_0.01_Prev,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                hubPar = "degree", hubQuant = 0.95,
                                normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(combined_analysis_0.05,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.1,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of Combined Datasets", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)

ps_EMP_genus <- tax_glom(ps_EMP, taxrank="Genus")
ps_EMP_1 <- core(ps_EMP_genus, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(ps_EMP_1) > 0
complete_data2 <- prune_samples(non_zero_samples, ps_EMP_1)
complete_data2 <- renameTaxa(complete_data2, 
                             pat = "<name>", 
                             substPat = "<subst_name>(<subst_R>)",
                             numDupli = "Genus")

adj_0.01_EMP <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP/degree_0.01_Prev.png")
net_0.01_EMP <- netConstruct(data=adj_0.01_EMP,
                              normMethod = "none", zeroMethod = "none",
                              sparsMethod = "none", dataType = "condDependence",
                              verbose = 3)
analysis_0.01_EMP <- netAnalyze(net_0.01_EMP,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                     hubPar = "degree", hubQuant = 0.95,
                                     normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(analysis_0.01_EMP,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.15,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of EMP Dataset", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)

Tara_genus <- tax_glom(Tara_phyloseq, taxrank ='Genus')
Tara_1 <- core(Tara_genus, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(Tara_1) > 0
complete_data2 <- prune_samples(non_zero_samples, Tara_1)
complete_data2 <- renameTaxa(complete_data2, 
                             pat = "<name>", 
                             substPat = "<subst_name>(<subst_R>)",
                             numDupli = "Genus")

adj_0.01_Tara <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Tara/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Tara/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Tara/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Tara/degree_0.01_Prev.png")
net_0.01_Tara <- netConstruct(data=adj_0.01_Prev,
                              normMethod = "none", zeroMethod = "none",
                              sparsMethod = "none", dataType = "condDependence",
                              verbose = 3)
analysis_0.01_Tara <- netAnalyze(net_0.01_Tara,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                     hubPar = "degree", hubQuant = 0.95,
                                     normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(analysis_0.01_Tara,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.15,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of Tara Oceans Dataset", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)

Qiita_2318_genus <- tax_glom(Qiita_2318_phyloseq, taxrank='Genus')
Qiita2318_1 <- core(Qiita_2318_genus, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(Qiita2318_1) > 0
complete_data2 <- prune_samples(non_zero_samples, Qiita2318_1)
complete_data2 <- renameTaxa(complete_data2, 
                             pat = "<name>", 
                             substPat = "<subst_name>(<subst_R>)",
                             numDupli = "Genus")

adj_0.01_2318 <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_2318/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_2318/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_2318/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_2318/degree_0.01_Prev.png")
net_0.01_2318 <- netConstruct(data=adj_0.01_2318,
                              normMethod = "none", zeroMethod = "none",
                              sparsMethod = "none", dataType = "condDependence",
                              verbose = 3)
analysis_0.01_2318 <- netAnalyze(net_0.01_2318,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                     hubPar = "degree", hubQuant = 0.95,
                                     normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(analysis_0.01_2318,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.20,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of Qiita 2318 Dataset", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)



Qiita_1552_genus <- tax_glom(Qiita_1552_phyloseq, taxrank='Genus')
Qiita1552_1 <- core(Qiita_1552_genus, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(Qiita1552_1) > 0
complete_data2 <- prune_samples(non_zero_samples, Qiita1552_1)
complete_data2 <- renameTaxa(complete_data2, 
                             pat = "<name>", 
                             substPat = "<subst_name>(<subst_R>)",
                             numDupli = "Genus")

adj_0.01_1552 <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_1552/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_1552/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_1552/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_1552/degree_0.01_Prev.png")
net_0.01_1552 <- netConstruct(data=adj_0.01_1552,
                              normMethod = "none", zeroMethod = "none",
                              sparsMethod = "none", dataType = "condDependence",
                              verbose = 3)
analysis_0.01_1552 <- netAnalyze(net_0.01_1552,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                 hubPar = "degree", hubQuant = 0.95,
                                 normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(analysis_0.01_1552,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.0,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of Qiita 1552 Dataset", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)


Atlantic_genus <- tax_glom(merge_atlantic, taxrank='Genus')
Atlantic_1 <- core(Atlantic_genus, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(Atlantic_1) > 0
complete_data2 <- prune_samples(non_zero_samples, Atlantic_1)
complete_data2 <- renameTaxa(complete_data2, 
                             pat = "<name>", 
                             substPat = "<subst_name>(<subst_R>)",
                             numDupli = "Genus")

adj_0.01_atlantic <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Atlantic_Ocean/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Atlantic_Ocean/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Atlantic_Ocean/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Atlantic_Ocean/degree_0.01_Prev.png")
net_0.01_atlantic <- netConstruct(data=adj_0.01_atlantic,
                              normMethod = "none", zeroMethod = "none",
                              sparsMethod = "none", dataType = "condDependence",
                              verbose = 3)
analysis_0.01_atlantic <- netAnalyze(net_0.01_atlantic,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                 hubPar = "degree", hubQuant = 0.95,
                                 normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(analysis_0.01_atlantic,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.0,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of Atlantic Ocean Dataset", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)


Qiita10178_genus <- tax_glom(Qiita_10178_phyloseq, taxrank='Genus')
Qiita10178_1 <- core(Qiita10178_genus, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(Qiita10178_1) > 0
complete_data2 <- prune_samples(non_zero_samples, Qiita10178_1)
complete_data2 <- renameTaxa(complete_data2, 
                             pat = "<name>", 
                             substPat = "<subst_name>(<subst_R>)",
                             numDupli = "Genus")

adj_0.01_10178 <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_10178/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_10178/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_10178/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/Qiita_10178/degree_0.01_Prev.png")
net_0.01_10178 <- netConstruct(data=adj_0.01_10178,
                                  normMethod = "none", zeroMethod = "none",
                                  sparsMethod = "none", dataType = "condDependence",
                                  verbose = 3)
analysis_0.01_10178 <- netAnalyze(net_0.01_10178,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                     hubPar = "degree", hubQuant = 0.95,
                                     normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(analysis_0.01_10178,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.15,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of Qiita 10178 Dataset", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)


EMP500_genus <- tax_glom(merge_EMP500, taxrank='Genus')
EMP500_1 <- core(EMP500_genus, detection=0, prevalence = 0.1, include.lowest = FALSE)
non_zero_samples <- sample_sums(EMP500_1) > 0
complete_data2 <- prune_samples(non_zero_samples, EMP500_1)
complete_data2 <- renameTaxa(complete_data2, 
                             pat = "<name>", 
                             substPat = "<subst_name>(<subst_R>)",
                             numDupli = "Genus")

adj_0.01_EMP500 <- node_gen(complete_data2, "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP500/node_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP500/edges_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP500/metrics_0.01_Prev.csv", "C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/EMP500/degree_0.01_Prev.png")
net_0.01_EMP500 <- netConstruct(data=adj_0.01_EMP500,
                               normMethod = "none", zeroMethod = "none",
                               sparsMethod = "none", dataType = "condDependence",
                               verbose = 3)
analysis_0.01_EMP500 <- netAnalyze(net_0.01_EMP500,centrLCC = FALSE,avDissIgnoreInf = TRUE,
                                  hubPar = "degree", hubQuant = 0.95,
                                  normDeg = TRUE, normBetw = TRUE, normEigen = TRUE, verbose=2)
#png("C:/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Results/SPIEC-EASI/Absolute_Abundance2/complete/network_image_1000.png")
plot_combined <- plot(analysis_0.01_EMP500,
                      sameLayout = TRUE,
                      # layout = "circle",
                      edgeFilter="threshold",
                      edgeFilterPar = 0.15,
                      layoutGroup = "union",
                      nodeColor = "cluster", 
                      rmSingles = "inboth",
                      nodeTransp = 0,
                      nodeSize = "degree",
                      highlightHubs = TRUE,
                      title1 = "Co-occurrence network of EMP500 Dataset", 
                      showTitle = TRUE,
                      cexNodes = 1,
                      cexHubs = 6,
                      cexHubLabels = 3,
                      cexTitle = 1.5,
                      mar =c(3, 4, 5, 2),
                      labelScale = TRUE,
                      labelFont = 2,
                      labels = NULL
)


