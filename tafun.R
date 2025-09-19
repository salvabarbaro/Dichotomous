library(IC2)
library(cluster)
library(dplyr)
library(parallel)

ta.fun <- function(df) {
  # Step 1: Gini decomposition
  df <- df %>% mutate(Approv = as.factor(Approv))  # Ensure grouping var is factor
  gini_decomp <- decompSGini(x = df$Rating, z = df$Approv)
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
#
## Example Ebert 2010
ex.ue <- data.frame(Candidate = LETTERS[1:4],
                   Rating = c(10, 20, 15, 15),
                   Approv = c(0, 0, 1, 1))

ta.fun(ex.ue)

ex.ue1 <- data.frame(Candidate = LETTERS[1:6],
                   Rating = c(.10, .23, .12, .15, .15, .15),
                   Approv = c(0, 0, 0, 1, 1, 1))
ta.fun(ex.ue1)

ex.ue2 <- data.frame(Candidate = LETTERS[1:6],
                   Rating = c(.14, .14, .17, .15, .15, .15),
                   Approv = c(0, 0, 0, 1, 1, 1))
ta.fun(ex.ue2)

due1.YL <- decompSGini(x = ex.ue1$Rating, z = as.factor(ex.ue1$Approv), decomp = "YL")
## G = 0.137, W = 0.096, B = 0.00
due1.BM <- decompSGini(x = ex.ue1$Rating, z = as.factor(ex.ue1$Approv), decomp = "BM")
## G = 0.137, W = 0.048, B = 00, overlap = 0.088
due1.A <- decompAtkinson(x = ex.ue1$Rating, z = as.factor(ex.ue1$Approv), epsilon = 1, decomp = "BDA")
## A = 0.0329, W = 0.0329, B = 0
due1.T <- decompGEI(x = ex.ue1$Rating, z = as.factor(ex.ue1$Approv))
## T = 0.034, W = 0.0344, B = 0
#######################################################################################
due2.YL <- decompSGini(x = ex.ue2$Rating, z = as.factor(ex.ue2$Approv), decomp = "YL")
## G = 0.033, W = 0.022, B = 0.00
due2.BM <- decompSGini(x = ex.ue2$Rating, z = as.factor(ex.ue2$Approv), decomp = "BM")
## G = 0.033, W = 0.011, B = 00, overlap = 0.002
due2.A <- decompAtkinson(x = ex.ue2$Rating, z = as.factor(ex.ue2$Approv), epsilon = 1, decomp = "BDA")
## A = 0.0021, W = 0.0021, B = 0
due2.T <- decompGEI(x = ex.ue2$Rating, z = as.factor(ex.ue2$Approv))
## T = 0.0021, W = 0.0021, B = 0
####################################################################################
## Check the following claim: if the preferences are threshold-dichotomous, then the BM-overlap is zero.
bmoverl.1 <- data.frame(Candidate = LETTERS[1:6],
                   Rating = c(.10, .16, .18, .2, .25, .29),
                   Approv = c(0, 0, 0, 1, 1, 1))
bm.1 <- decompSGini(x = bmoverl.1$Rating, z = as.factor(bmoverl.1$Approv), decomp = "BM")
# G = 0.175, W = 0.048, B = 0.127, ### W+B = G
bmoverl.2 <- data.frame(Candidate = LETTERS[1:6],
                   Rating = c(.10, .16, .2, .18, .25, .29),
                   Approv = c(0, 0, 0, 1, 1, 1))
bm.2 <- decompSGini(x = bmoverl.2$Rating, z = as.factor(bmoverl.2$Approv), decomp = "BM")
# G = 0.175, W = 0.059, B = 0.11, Overlap = 0.0056, ### W+B = G
bm.2$decomp$within + bm.2$decomp$between + bm.2$decomp$overlap
bm.2$ineq
#
yl.2 <- decompSGini(x = bmoverl.2$Rating, z = as.factor(bmoverl.2$Approv), decomp = "YL")
# G = 0.175, W = 0.048, B = 0.127, ### W+B = G
as.numeric(yl.2$ineq)[1] ## Gini
yl.2$decomp$within
yl.2$decomp$between
yl.2$decomp$stratif
yl.2$decomp$within + yl.2$decomp$between + yl.2$decomp$stratif


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
