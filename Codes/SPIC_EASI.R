############################################################
## SPIEC-EASI Network Analysis Pipeline
############################################################

library(SpiecEasi)
library(igraph)
library(dplyr)

############################################################
## INPUT
############################################################
combined_metadata_16s_Control <- combined_metadata_16s[combined_metadata_16s$study_condition == "Control",]
combined_metadata_16s_Diseased <- combined_metadata_16s[combined_metadata_16s$study_condition == "Diseased",]

# Replace with your dataset
count_matrix_control <- species[rownames(species) %in% combined_metadata_16s_Control$sample_id,]
count_matrix_diseased <- species[rownames(species) %in% combined_metadata_16s_Diseased$sample_id,]

############################################################
## PREPROCESSING
############################################################

# Replace NA values with zero
count_matrix_control[is.na(count_matrix_control)] <- 0
count_matrix_diseased[is.na(count_matrix_diseased)] <- 0

# Remove samples with zero counts
count_matrix_control <- as.matrix(count_matrix_control[rowSums(count_matrix_control) > 0, ])
count_matrix_diseased <- as.matrix(count_matrix_diseased[rowSums(count_matrix_diseased) > 0, ])


############################################################
## RUN SPIEC-EASI for Control
############################################################

se.mb_control <- spiec.easi(
  count_matrix_control,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.params = list(
    rep.num = 50,
    seed = 100
  )
)


############################################################
## RUN SPIEC-EASI for diseased
############################################################

se.mb_diseased <- spiec.easi(
  count_matrix_diseased,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.params = list(
    rep.num = 50,
    seed = 100
  )
)


############################################################
## CREATE NETWORK for control
############################################################

adj_control <- getRefit(se.mb_control)

rownames(adj_control) <- colnames(count_matrix_control)
colnames(adj_control) <- colnames(count_matrix_control)

g_control <- graph_from_adjacency_matrix(
  adj_control,
  mode = "undirected",
  diag = FALSE
)


############################################################
## CREATE NETWORK for Diseased
############################################################

adj_diseased <- getRefit(se.mb_diseased)

rownames(adj_diseased) <- colnames(count_matrix_diseased)
colnames(adj_diseased) <- colnames(count_matrix_diseased)

g_diseased <- graph_from_adjacency_matrix(
  adj_diseased,
  mode = "undirected",
  diag = FALSE
)

############################################################
## NETWORK STATISTICS
############################################################

cat("\n----------------------------\n")
cat("NETWORK SUMMARY\n")
cat("----------------------------\n")

cat("Nodes:", vcount(g_control), "\n")
cat("Edges:", ecount(g_control), "\n")
cat("Density:", edge_density(g_control), "\n")
cat("Clustering coefficient:", transitivity(g_control), "\n")
cat("Diameter:", diameter(g_control), "\n")

comp_control <- components(g_control)

cat("Connected components:", comp_control$no, "\n")
cat("Largest component:", max(comp_control$csize), "nodes\n")


############################################################

cat("\n----------------------------\n")
cat("NETWORK SUMMARY\n")
cat("----------------------------\n")

cat("Nodes:", vcount(g_diseased), "\n")
cat("Edges:", ecount(g_diseased), "\n")
cat("Density:", edge_density(g_diseased), "\n")
cat("Clustering coefficient:", transitivity(g_diseased), "\n")
cat("Diameter:", diameter(g_diseased), "\n")

comp_diseased <- components(g_diseased)

cat("Connected components:", comp_diseased$no, "\n")
cat("Largest component:", max(comp_diseased$csize), "nodes\n")

############################################################
## COMMUNITY DETECTION for control
############################################################

comm_control <- cluster_louvain(g_control)

cat("\nCommunities:", length(sizes(comm_control)), "\n")
print(sizes(comm_control))

V(g_control)$Community <- membership(comm_control)

############################################################
## CENTRALITY ANALYSIS for control
############################################################

hub_df_control <- data.frame(
  Taxon = V(g_control)$name,
  Degree = degree(g_control),
  Betweenness = betweenness(g_control),
  Closeness = closeness(g_control),
  Eigenvector = eigen_centrality(g_control)$vector,
  Community = membership(comm_control)
)

hub_df_control <- hub_df_control %>%
  arrange(desc(Degree))

cat("\nTop 20 Hub Taxa\n")
print(head(hub_df_control,20))


############################################################
## COMMUNITY DETECTION for diseased
############################################################

comm_diseased <- cluster_louvain(g_diseased)

cat("\nCommunities:", length(sizes(comm_diseased)), "\n")
print(sizes(comm_diseased))

V(g_diseased)$Community <- membership(comm_diseased)

############################################################
## CENTRALITY ANALYSIS for diseased
############################################################

hub_df_diseased <- data.frame(
  Taxon = V(g_diseased)$name,
  Degree = degree(g_diseased),
  Betweenness = betweenness(g_diseased),
  Closeness = closeness(g_diseased),
  Eigenvector = eigen_centrality(g_diseased)$vector,
  Community = membership(comm_diseased)
)

hub_df_diseased <- hub_df_diseased %>%
  arrange(desc(Degree))

cat("\nTop 20 Hub Taxa\n")
print(head(hub_df_diseased,20))




############################################################
## TOP HUBS PER COMMUNITY
############################################################

community_hubs <- hub_df %>%
  group_by(Community) %>%
  slice_max(Degree,
            n = 5,
            with_ties = FALSE)

print(community_hubs)

############################################################
## EXTRACT LARGEST CONNECTED COMPONENT for control
############################################################

largest <- induced_subgraph(
  g_,
  which(comp$membership == which.max(comp$csize))
)

cat("\nLargest Component\n")
cat("Nodes:",vcount(largest),"\n")
cat("Edges:",ecount(largest),"\n")

############################################################
## EXPORT TO CYTOSCAPE
############################################################

write_graph(
  g,
  "SPIEC_network.graphml",
  format="graphml"
)

write_graph(
  largest,
  "SPIEC_network_largest_component.graphml",
  format="graphml"
)

############################################################
## EXPORT NODE ATTRIBUTES
############################################################

write.csv(
  hub_df,
  "Node_attributes_diseased.csv",
  row.names=FALSE
)

############################################################
## EXPORT COMMUNITY MEMBERSHIP
############################################################

community_df <- data.frame(
  Taxon = V(g)$name,
  Community = membership(comm)
)

write.csv(
  community_df,
  "Community_assignment.csv",
  row.names=FALSE
)

############################################################
## EXPORT EDGE LIST
############################################################

edge_df <- as_data_frame(g, what="edges")

write.csv(
  edge_df,
  "Edge_list.csv",
  row.names=FALSE
)

############################################################
## POSITIVE / NEGATIVE ASSOCIATIONS
############################################################

beta <- symBeta(getOptBeta(se.mb), mode="maxabs")

positive_edges <- sum(beta > 0)/2
negative_edges <- sum(beta < 0)/2

cat("\nPositive edges:",positive_edges,"\n")
cat("Negative edges:",negative_edges,"\n")

############################################################
## FINISHED
############################################################
