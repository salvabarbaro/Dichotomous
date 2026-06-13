library(dplyr)
library(cluster)

x <- c(1.3, 1.9, 1.5, 1.9, 1.6, 1.9)

df <- tibble(
  id = seq_along(x),
  x = x
)

set.seed(123)

km <- kmeans(df$x, centers = 2, nstart = 25)

df <- df %>%
  mutate(cluster = km$cluster)

# Distance matrix
D <- as.matrix(dist(df$x))

# Cohesion a(i): average distance to points in own cluster
# Separation b(i): smallest average distance to another cluster
sil_info <- df %>%
  rowwise() %>%
  mutate(
    cohesion = {
      own <- which(df$cluster == cluster & df$id != id)
      if (length(own) == 0) 0 else mean(D[id, own])
    },
    separation = {
      other_clusters <- setdiff(unique(df$cluster), cluster)
      min(sapply(other_clusters, function(cl) {
        mean(D[id, df$cluster == cl])
      }))
    },
    silhouette = ifelse(
      max(cohesion, separation) == 0,
      0,
      (separation - cohesion) / max(cohesion, separation)
    )
  ) %>%
  ungroup()

sil_info
stargazer::stargazer(sil_info, summary = F)


sil <- silhouette(km$cluster, dist(df$x))
sil

mean(sil[, "sil_width"])
plot(sil)


###########################################################################
D <- dist(x)

avg_sil <- function(k, x, D) {
  
  km <- kmeans(
    matrix(x, ncol = 1),
    centers = k,
    nstart = 50
  )
  
  sil <- silhouette(km$cluster, D)
  
  mean(sil[, "sil_width"])
}

results <- data.frame(
  k = 2:4,
  avg_silhouette = sapply(2:4, avg_sil, x = x, D = D)
)

results


library(ggplot2)

ggplot(results, aes(x = k, y = avg_silhouette)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(breaks = 1:4) +
  labs(
    x = "Number of clusters (k)",
    y = "Average silhouette width"
  ) +
  theme_bw(base_size = 18)


fviz_nbclust(r, kmeans, method = "silhouette", k.max = 4)['data'] 
