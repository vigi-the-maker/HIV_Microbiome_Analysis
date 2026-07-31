# Set thresholds
min_prevalence <- 0.30
min_total_count <- 400

# Calculate prevalence and total abundance
prevalence <- colSums(combined_species_profile > 0)
total_counts <- colSums(combined_species_profile)

# Keep taxa meeting both criteria
keep <- prevalence >= (min_prevalence * nrow(combined_species_profile)) &
  total_counts >= min_total_count

# Filter the count matrix
combined_species_profile_filtered <- combined_species_profile[,keep]

# Check dimensions
dim(count_matrix)
dim(count_matrix_filtered)