library(IC2)
library(cluster)
library(dplyr)
library(parallel)

ta.fun <- function(df) {
  # Step 1: Gini decomposition
  df <- df %>% mutate(Approv = as.factor(Approv))  # Ensure grouping var is factor
  gini_decomp <- decompSGini(x = df$Rating, z = df$Approv, decomp = "YL")
  gini_within <-  gini_decomp$decomp$within
  gini_between <- gini_decomp$decomp$between
  WDP <- ifelse(gini_within < gini_between, "Yes","No")

  # Step 2: Optimal number of clusters
  X <- df %>% select(Approv, Rating) %>%
    mutate(Approv = as.numeric(as.character(Approv)))

  # Get number of unique observations
  num_unique <- nrow(unique(X))
  max_k <- min(6, num_unique - 1)

  if (max_k >= 2) {
    dist_mat <- dist(X)
    silhouette_scores <- sapply(2:max_k, function(k) {
      km_res <- kmeans(X, centers = k, nstart = 25)
      sil <- silhouette(km_res$cluster, dist_mat)
      mean(sil[, 3])
    })
    optimal_k <- which.max(silhouette_scores) + 1  # shift because k starts at 2
  } else {
    optimal_k <- NA  # Not enough variation for clustering
  }

  # Step 3: Approval separation test
  approved_min <- min(df$Rating[df$Approv == 1], na.rm = TRUE)
  disapp_max   <- max(df$Rating[df$Approv == 0], na.rm = TRUE)
  threshold <- ifelse(disapp_max <= approved_min, "Yes", "No")
  approval_test <- as.integer(disapp_max <= approved_min)

  # Return result
  result <- data.frame(
    Gini_Within = gini_within,
    Gini_Between = gini_between,
    WDP = WDP,
    Optimal_k = optimal_k,
    Approval_Sep = threshold
  )

  return(result)
}




#ex1 <- data.frame(Candidate = LETTERS[1:6],
#                  Approv = c(0, 0, 0, 1, 1, 1),
#                  Rating = c(1.3, 1.9, 1.5, 1.6, 1.9, 1.6))
ex1 <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 0, 1, 1, 1),
                  Rating = c(1.3, 1.3, 1.3, 1.8, 1.8, 1.8))
ex2 <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 0, 1, 1, 1),
                  Rating = c(1.2, 1.3, 1.4, 1.7, 1.8, 1.9))
ex3 <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 0, 1, 1, 1),
                  Rating = c(1.3, 1.5, 1.6, 1.7, 1.7, 1.9))
ex4 <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 0, 1, 1, 1),
                  Rating = c(1.3, 1.9, 1.5, 1.6, 1.9, 1.6))
ex5 <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 1, 0, 1, 1),
                  Rating = c(1.1, 1.1, 1.4, 1.5, 1.9, 1.9))
## Note: with funny ex5 yields three clusters.

ta.fun(ex5)

## Example from the theoretical section to explain the bleeding persuasiveness of the threshold approach
ex.t <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 0, 0, 0, 1),
                  Rating = c(.1, .2, .4, .6, .8, 0.9))
ta.fun(ex.t)
##########################################################################################################

ex.list <- list(ex1 = ex1, ex2 = ex2, ex3 = ex3, ex4 = ex4, ex5 = ex5)
tb1 <- lapply(ex.list, ta.fun) %>% bind_rows(., .id = "Case")
tb1


#########################################################################
### Generate 10,000 data frames like exX
set.seed(55234)

generate.fun <- function(m) {
  data.frame(
    Candidate = LETTERS[1:m],
    Approv = sample(c(0, 1), size = m, replace = TRUE),
    Rating = runif(m, min = 1, max = 2)  # Uniform real numbers between 1 and 2
  )
}

# Generate 10,000 data frames in a list
df.list <- replicate(10000, generate.fun(m=6), simplify = FALSE)

# Optional: name the list elements
names(df.list) <- paste0("ex", seq_along(df.list))

tb.repl <- mclapply(df.list, ta.fun, mc.cores = 8) %>% bind_rows(., .id = "Case")
tb.repl
k3 <- tb.repl %>% filter(., Optimal_k == 3)
k3 %>% filter(., Approval_Sep == "Yes" & WDP == "No")
df.list$ex9479


k3<- tb.repl %>% filter(Optimal_k == 3)
k3.filtered <- k3 %>% filter(Approval_Sep == "Yes", WDP == "No")

# Step 2: Get the case names
case_ids <- k3.filtered$Case

# Step 3: Check Approv count in original data frames
case_ids_filtered <- case_ids[sapply(case_ids, function(id) {
  df <- df.list[[id]]
  num_approved <- sum(df$Approv == 1)
  return(num_approved >= 2 && num_approved <= 4)
})]

# Step 4: Output
print(case_ids_filtered)
df.list$ex7487






########## aux

gini_decomp <- decompSGini(x = ex1$Rating, z = factor(ex1$Approv), decomp = "YL")
  gini_within <-  gini_decomp$decomp$within
  gini_between <- gini_decomp$decomp$between
## gini: 0.0806
gini_within + gini_between
gini_within
gini_between


############################
#Theil-T
exB <- data.frame(Candidate = LETTERS[1:6],
                  Approv = c(0, 0, 0, 1, 1, 1),
                  Rating = c(130, 140, 135, 170, 180, 190))


theil.fun <- function(df){
  df <- df %>% mutate(Approv = as.factor(Approv))
  dec <- decompGEI(x = df$Rating, z = df$Approv, alpha = 1, ELMO = FALSE)
  th.total <- as.numeric(dec$ineq$index)
  th.within <- as.numeric(dec$decomp$within)
  th.between <- as.numeric(dec$decomp$between)
  return(c(th.total, th.within, th.between))
}
theil.fun(exB)

mean(exB$Rating)

helpfun<- function(v){
  ka = v/157.5 * log(v/157.5)
  return(ka)
}


decB <- decompGEI(x = exB$Rating, z = as.factor(exB$Approv), alpha = 1, ELMO = F)


#############################
## Atkinson-Index, epsilon = 1
atk_decomp <- function(df, value_col = "Rating", group_col = "Approv") {
  x <- df[[value_col]]
  g <- df[[group_col]]
  m <- length(x)
  
  mu <- mean(x)
  gmean <- exp(mean(log(x)))
  atk_total <- 1 - gmean / mu
  
  # Group-level stats
  groups <- split(x, g)
  mu_g <- sapply(groups, mean)
  gmean_g <- sapply(groups, function(v) exp(mean(log(v))))
  m_g <- sapply(groups, length)
  
  # weights
  s_g <- (m_g * mu_g) / (m * mu)
  
  # subgroup Atkinson
  atk_g <- 1 - gmean_g / mu_g
  
  # within and between
  atk_within <- sum(s_g * atk_g)
  atk_between <- sum(s_g * (gmean_g / mu_g)) - (gmean / mu)
  
  list(
    total   = atk_total,
    within  = atk_within,
    between = atk_between,
    check   = atk_within + atk_between
  )
}
atk_decomp(df = exB)


atk <- decompAtkinson(x = exB$Rating, z = as.factor(exB$Approv), decomp = "BDA", epsilon = 1) 
as.numeric(atk$ineq$index) # total
as.numeric(atk$decomp$within) # within
as.numeric(atk$decomp$between) # between
