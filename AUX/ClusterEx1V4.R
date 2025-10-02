library(cluster)
library(dplyr)

## kmeans
ex1 <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 0, 1, 1, 1),
                  Rating = c(1.3, 1.9, 1.5, 1.6, 1.9, 1.6))

# Use only numeric features
X <- ex1 %>% select(Approv, Rating)

# Compute distance matrix (Euclidean)
dist_mat <- dist(X)

# Initialize results
silhouette_results <- data.frame(k = integer(), avg_sil_width = numeric())

# Loop over k = 1 to 3
for (k in 1:3) {
  if (k == 1) {
    # Silhouette is undefined for k = 1
    silhouette_results <- rbind(silhouette_results, data.frame(k = 1, avg_sil_width = NA))
  } else {
    # Run k-means clustering
    kmeans_res <- kmeans(X, centers = k, nstart = 25)
    
    # Compute silhouette
    sil <- silhouette(kmeans_res$cluster, dist_mat)
    
    # Compute average silhouette width
    avg_width <- mean(sil[, 3])
    
    # Store results
    silhouette_results <- rbind(silhouette_results, data.frame(k = k, avg_sil_width = avg_width))
  }
}

# Print results
print(silhouette_results)

