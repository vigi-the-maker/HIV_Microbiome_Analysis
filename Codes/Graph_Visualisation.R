############################################################
## VISUALIZATION USING QGRAPH
############################################################

library(qgraph)
library(RColorBrewer)
library(scales)



############################################################
## BASIC NETWORK for control
############################################################
# Community membership
membership_control <- membership(comm_control)

# Size of each community
community_sizes_control <- sizes(comm_control)

# Communities with more than one node
keep_communities_control <- names(community_sizes_control[community_sizes_control > 1])

# Vertices to keep
keep_vertices_control <- V(g_control)[membership_control %in% keep_communities_control]

# Subgraph
g_filtered_control <- induced_subgraph(g_control, keep_vertices_control)
membership_filtered_control <- membership(comm_control)[V(g_filtered_control)$name]
groups_control <- as.factor(membership_filtered_control)
deg_control <- degree(g_filtered_control)
community_sizes_control_filtered <- sort(community_sizes_control[community_sizes_control > 1],decreasing = TRUE)
top10 <- as.numeric(names(community_sizes_control_filtered)[1:10])

ec <- eigen_centrality(g_filtered_control)$vector

rep_species <- sapply(top10, function(comm) {
  nodes <- names(membership_filtered_control[membership_filtered_control == comm])
  nodes[which.max(deg_control[nodes])]
  
})



hub_nodes_control <- names(sort(deg, decreasing = TRUE))[1:10]


labels <- rep("", vcount(g_filtered_control))
names(labels) <- V(g_filtered_control)$name

pretty_labels <- gsub("_", " ", rep_species)
pretty_labels <- paste(pretty_labels, "Hub")

labels[rep_species] <- pretty_labels



node_colours_control <- rainbow(length(levels(groups_control)))[groups_control]
qgraph(
  as_adjacency_matrix(g_filtered_control),
  labels = labels,
  color = node_colours,
  vsize = scales::rescale(deg_control, to = c(4, 20)),
  layout = "spring",
  label.cex = 1.2,
  label.font = 2,
  title="Control Gut microbiome Network"
)



############################################################
## BASIC NETWORK for DISEASED
############################################################
# Community membership
membership_diseased <- membership(comm_diseased)

# Size of each community
community_sizes_diseased <- sizes(comm_diseased)

# Communities with more than one node
keep_communities_diseased <- names(community_sizes_diseased[community_sizes_diseased > 1])

# Vertices to keep
keep_vertices_diseased <- V(g_diseased)[membership_diseased %in% keep_communities_diseased]

# Subgraph
g_filtered_diseased <- induced_subgraph(g_diseased, keep_vertices_diseased)
membership_filtered_diseased <- membership(comm_diseased)[V(g_filtered_diseased)$name]
groups_diseased <- as.factor(membership_filtered_diseased)
deg_diseased <- degree(g_filtered_diseased)
community_sizes_diseased_filtered <- sort(community_sizes_diseased[community_sizes_diseased > 1],decreasing = TRUE)
top10_diseased <- as.numeric(names(community_sizes_diseased_filtered)[1:10])
ec_diseased <- eigen_centrality(g_filtered_diseased)$vector
community_sizes_filtered <- sort(
  table(membership_filtered_diseased),
  decreasing = TRUE
)

top10_diseased <- as.numeric(names(community_sizes_filtered)[1:10])
rep_species_diseased <- sapply(top10_diseased, function(comm) {
  nodes <- names(membership_filtered_control[membership_filtered_diseased == comm])
  nodes[which.max(deg_diseased[nodes])]
  
})



hub_nodes_control <- names(sort(deg, decreasing = TRUE))[1:10]


labels_diseased <- rep("", vcount(g_filtered_diseased))
names(labels_diseased) <- V(g_filtered_diseased)$name

pretty_labels <- gsub("_", " ", rep_species_diseased)
pretty_labels <- paste(pretty_labels, "Hub")

labels_diseased[rep_species_diseased] <- pretty_labels



node_colours_diseased <- rainbow(length(levels(groups_diseased)))[groups_diseased]
qgraph(
  as_adjacency_matrix(g_filtered_diseased),
  labels = labels_diseased,
  color = node_colours_diseased,
  vsize = scales::rescale(deg_diseased, to = c(4, 20)),
  layout = "spring",
  label.cex = 1.2,
  label.font = 2,
  title = "HIV Gut Microbiome Association Network"

)

############################################################
## SAVE HIGH-RESOLUTION PNG
############################################################

png(
  "SPIEC_network_qgraph.png",
  width = 5000,
  height = 5000,
  res = 600
)

qgraph(
  adj,
  layout = "spring",
  labels = labels,
  color = node_colours,
  vsize = rescale(deg,to=c(4,12)),
  edge.color = "grey70",
  label.cex = 0.8,
  border.color = "black"
)

dev.off()

############################################################
## SAVE PDF
############################################################

pdf(
  "SPIEC_network_qgraph.pdf",
  width = 14,
  height = 14
)

qgraph(
  adj,
  layout = "spring",
  labels = labels,
  color = node_colours,
  vsize = rescale(deg,to=c(4,12)),
  edge.color = "grey70",
  label.cex = 0.8,
  border.color = "black"
)

dev.off()
