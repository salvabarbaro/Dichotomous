### Graz Dataset
### I. Read Data and prepare data
setwd("~/Documents/Research/Dichotomous")
#########################################
library(dplyr)
library(cluster)
library(factoextra)
library(ggplot2)
library(ineq)
library(rstatix)
library(tidyr)
library(IC2)
library(parallel)
library(gtsummary)
library(readxl)
#########################################
# Section 0: Read data and data handling
austria.df <- read_excel("Data/Steirische_LTW_2019_Daten_Barbaro.xlsx") %>%
  rename(id = lfdn) %>%
  mutate(Gender = factor(Qb, 
                         levels = c(1, 2, 3), 
                         labels = c("male", "female", "diverse")),
         Gender = replace(Gender, Qb %in% c(98, 99), NA)) %>%
  mutate(Age = factor(Qc,
                      levels = 1:7,
                      labels = c("16to19", "20to29", "30to39", "40to49",
                                 "50to59", "60to69", "70+")),
         Age = replace(Age, Qc %in% c(98, 99), NA)) %>%
  mutate(Educ = factor(Qd,
                       levels = 1:6,
                       labels = c("Pflichtschule", "Pflichtschule+Lehre",
                                  "Mittlere Schule", "AHS/BHS", "Matura", "HigherEduc") ),
         Educ = replace(Educ, Qd %in% c(98,99), NA) ) %>%
  mutate(across(starts_with("Q1"), ~ replace(., . == 99, NA))) %>%
  rename(Approv.SPÖ = Q1_1,
         Approv.ÖVP = Q1_2,
         Approv.FPÖ = Q1_3,
         Approv.Green = Q1_4,
         Approv.KPÖ = Q1_5,
         Approv.NEOS = Q1_6) %>%
  mutate(across(starts_with("Q5"), ~ replace(., . %in% 97:99, NA))) %>%
  rename(Rat.SPÖ = Q5_1,
         Rat.ÖVP = Q5_2,
         Rat.FPÖ = Q5_3,
         Rat.Green = Q5_4,
         Rat.KPÖ = Q5_5,
         Rat.NEOS = Q5_6) %>%
  mutate(across(starts_with("Q6"), ~ replace(., . %in% 97:99, NA))) %>%
  rename(Trich.SPÖ = Q6_1,
         Trich.ÖVP = Q6_2,
         Trich.FPÖ = Q6_3,
         Trich.Green = Q6_4,
         Trich.KPÖ = Q6_5,
         Trich.NEOS = Q6_6)

# Data frame for Theil analysis
austria_long.df <- austria.df %>%
  select(id, starts_with("Approv."), starts_with("Rat."))  %>%
  pivot_longer(cols = starts_with("Rat."), names_to = "Party", values_to = "Rating") %>%
  pivot_longer(cols = starts_with("Approv."), names_to = "Approval_Party", values_to = "Approval") %>%
  filter(substring(Party, 5) == substring(Approval_Party, 8)) %>%  # Ensure party names match
  mutate(Party = gsub("^Rat\\.", "", Party)) %>%
  select(., -c("Approval_Party"))

# Function to transform rating values on a [2,3]-scale
#Transformation on a [2,3]-scale - define function
scale_to_range <- function(x, new_min, new_max) {
  old_min <- min(x, na.rm = T)
  old_max <- max(x, na.rm = T)
  scaled_x <- ((x - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
  return(scaled_x)
}
# Transformation... - apply function
newratings <- scale_to_range(x = austria_long.df$Rating, 2,3)
# Transformation ... - replace original values with new rating values
austria_theil.df <- austria_long.df %>% 
  mutate(Rating = newratings)

# Data frame for cluster analysis (na_omit of Rating)
austria_cluster.df <- austria_theil.df %>% 
  group_by(id) %>%
  filter(all(!is.na(Rating))) %>%
  ungroup()

#############################################################
### Section 3: Decomposition analysis
ic2decomp.fun <- function(i) {
  tryCatch({
    df <- austria_theil.df %>% 
      filter(id %in% i) %>% 
      mutate(Approval = factor(Approval)) 
    
    dec <- decompGEI(x = df$Rating, z = df$Approval, alpha = 1, ELMO = FALSE)
    th.total <- as.numeric(dec$ineq$index)
    th.within <- as.numeric(dec$decomp$within)
    th.between <- as.numeric(dec$decomp$between)
    # Compute Gini decomposition
    gin <- decompSGini(x = df$Rating, z = df$Approval, decomp = "YL", ELMO = FALSE)
    gi.total <- as.numeric(gin$ineq$index)
    gi.within <- as.numeric(gin$decomp$within)
    gi.between <- as.numeric(gin$decomp$between)
    # Compute Atkinson decomposition
    atk <- decompAtkinson(x = df$Rating, z = df$Approval, decomp = "BDA", epsilon = 1)
    atk.total <- as.numeric(atk$ineq$index)
    atk.within <- as.numeric(atk$decomp$within)
    atk.between <- as.numeric(atk$decomp$between)
    # Store results in a data frame
    res.df <- data.frame(
      id = unique(df$id),
      Theil.tot = th.total,
      Theil.within = th.within,
      Theil.between = th.between,
      Gini.tot = gi.total,
      Gini.within = gi.within,
      Gini.between = gi.between,
      Atkinson.tot = atk.total,
      Atkinson.within = atk.within,
      Atkinson.between = atk.between
    )
    return(res.df)
  }, error = function(e) {
    message(paste("Skipping id:", i, "due to error:", e$message))
  })
}

# Running the function using lapply, removing NULL results
ic2res <- lapply(unique(austria_theil.df$id), ic2decomp.fun)
ic2res.df <- ic2res %>% do.call(rbind, .) %>%
  mutate(WDP.Theil = ifelse(Theil.within < Theil.between, 1, 0),
         WDP.Gini =  ifelse(Gini.within < Gini.between, 1, 0),
         WDP.Atkinson = ifelse(Atkinson.within < Atkinson.between, 1, 0))
head(ic2res.df)

compute_wdp_shares <- function(df) {
  wdp_share <- function(column) {
    valid_values <- column[!is.na(column)]  # Remove NA values
    if (length(valid_values) == 0) return(NA)  # Avoid division by zero
    return(sum(valid_values) / length(valid_values))
  }
  
  return(c(
    Theil = wdp_share(df$WDP.Theil),
    Gini = wdp_share(df$WDP.Gini),
    Atkinson = wdp_share(df$WDP.Atkinson)
  ))
}
## Values for Table 2:
compute_wdp_shares(ic2res.df)

####################################################################
## Bootstrap
# Function to compute WDP shares
compute_wdp_shares <- function(df) {
  wdp_share <- function(column) {
    valid_values <- column[!is.na(column)]  # Remove NA values
    if (length(valid_values) == 0) return(NA)  # Avoid division by zero
    return(sum(valid_values) / length(valid_values))
  }
  
  return(c(
    Theil = wdp_share(df$WDP.Theil),
    Gini = wdp_share(df$WDP.Gini),
    Atkinson = wdp_share(df$WDP.Atkinson)
  ))
}

# Custom bootstrapping function
cluster_bootstrap <- function(df, id_col, R = 1000) {
  unique_ids <- unique(df[[id_col]])  # Unique id values
  boot_results <- matrix(NA, nrow = R, ncol = 3)  # Store bootstrap results
  
  for (r in 1:R) {
    sampled_ids <- sample(unique_ids, replace = TRUE)  # Resample IDs
    boot_df <- df %>% filter(id %in% sampled_ids)  # Extract corresponding rows
    
    boot_results[r, ] <- compute_wdp_shares(boot_df)  # Compute WDP shares
  }
  
  colnames(boot_results) <- c("Theil", "Gini", "Atkinson")
  return(as.data.frame(boot_results))
}

# Running the clustered bootstrap
set.seed(55234)  
boot_results <- cluster_bootstrap(ic2res.df, id_col = "id", R = 1000)
# Compute confidence intervals by percentile method
apply(boot_results, 2, quantile, probs = c(0.025, 0.975))  # 95% CI

##################################################################
#### SECTION 4 Cluster-analytical assessment
## Cluster analysis
## present: austria_cluster.df
# Step 1: We remove all id's without variance in the Rating
# Function working.ids selects the appropriate id-values 
working.ids <- austria_cluster.df %>% 
  group_by(id) %>% 
  reframe(l = var(Rating, na.rm = T)) %>%
  filter(., l > 0) 
kmax.value <- 5  # for Graz, 10 for Grenoble

fv.fun <- function(i, df){
  ro <- df %>% filter(., id %in% i)
  epsilons <- sample(seq(-.1, .1, .001), 6, replace = F)
  rat.temp <- ro$Rating + epsilons
  sc.rat <- scale(rat.temp)
  silh.values <- fviz_nbclust(sc.rat, kmeans, 
                              method = "silhouette", 
                              k.max = kmax.value)[["data"]][["y"]]
  df.max <- data.frame(nbcluster = 1:(kmax.value),
                       silh.scores = silh.values)
  optk <- df.max[,1][which.max(df.max[,2])]
  return(optk)
}

optclust.list <- mclapply(working.ids$id, fv.fun, 
                          df = austria_cluster.df, mc.cores = 14)
optclust.df <- data.frame(id = working.ids$id, 
                          optk = unlist(optclust.list))
## Data for the Table in Section 4
table(optclust.df$optk) / nrow(optclust.df)


## Second Approach: Select the id's with kopt == 2 and check the silhouette scores
optclust2.df <- optclust.df %>% filter(., optk == 2)
sil.fun <- function(i, df){
  ro <- df %>% filter(., id %in% i)
  epsilons <- sample(seq(-.1, .1, .001), 6, replace = F)
  rat.temp <- ro$Rating + epsilons
  sc.rat <- scale(rat.temp)
  #  silh.values <- fviz_nbclust(sc.rat, kmeans, 
  #                              method = "silhouette", 
  #                              k.max = 2)[["data"]][["y"]]
  km.res <- kmeans(sc.rat, centers = 2)
  silh.values <- silhouette(km.res$cluster, dist = dist(sc.rat)) %>%
    as.data.frame()
  mn.silh = mean(silh.values$sil_width)
  return(mn.silh)
}

k2.ids <- working.ids %>% filter(., id %in% optclust2.df$id)
sil.values <- lapply(k2.ids$id, sil.fun, df = austria_cluster.df)
silh.df <- data.frame(id = optclust2.df$id,
                      silh.values = unlist(sil.values)) %>%
  mutate(
    silhouette_category = case_when(
      silh.values > 0.70 ~ "strong",
      silh.values > 0.50 ~ "moderate",
      silh.values > 0.25 ~ "weak",
      silh.values > 0    ~ "poor",
      silh.values <= 0   ~ "incorrect",
      is.na(silh.values) ~ NA_character_  # Keep NA if silhouette is missing
    )
  )
table(silh.df$silhouette_category)
table(silh.df$silhouette_category) / nrow(silh.df)

### BEGIN FANNY #############################################################
### ROBUSTNESS CHECK VIA FANNY
# First step: find out the optimal cluster number for each respondent
# Hint: because n=6 --> kmax.value = 2, this analysis is senseless for the Graz
# dataset. However, to ensure compatability to the other replication files, 
# we run the analysis nevertheless.
kmax.value = 2 # for Graz, higher in the others
faopt.fun <- function(i, df){
  ro <- df %>% filter(., id %in% i)
  epsilons <- sample(seq(-.1, .1, .001), 6, replace = F)
  rat.temp <- ro$Rating + epsilons
  sc.rat <- scale(rat.temp)
  silh.values <- fviz_nbclust(sc.rat, fanny, 
                              method = "silhouette", 
                              k.max = kmax.value)[["data"]][["y"]]
  df.max <- data.frame(nbcluster = 1:(kmax.value),
                       silh.scores = silh.values)
  optk <- df.max[,1][which.max(df.max[,2])]
  return(optk)
}
optclust_fa.list <- mclapply(working.ids$id, faopt.fun, 
                             df = austria_cluster.df, mc.cores = 14)
optclust_fa.df <- data.frame(id = working.ids$id, 
                             optk = unlist(optclust_fa.list))
## Second step: Evaluate the silhouette width for each respondent. 
## Values in Table!
fasil.fun <- function(i){
  df <- austria_cluster.df %>% filter(., id %in% i)
  fan = fanny(scale(df$Rating), k = 2)
  silh = fan$silinfo$avg.width
  return(silh)
}

fa.sil.list <- lapply(unique(working.ids$id), FUN=fasil.fun)
fa.sil.df <- data.frame(id = working.ids$id,
                        silh.values = unlist(fa.sil.list)) %>%
  mutate(silhouette_category = case_when(
    silh.values > 0.70 ~ "strong",
    silh.values > 0.50 ~ "moderate",
    silh.values > 0.25 ~ "weak",
    silh.values > 0    ~ "poor",
    silh.values <= 0   ~ "incorrect",
    is.na(silh.values) ~ NA_character_  # Keep NA if silhouette is missing
  ))
table(fa.sil.df$silhouette_category)
table(fa.sil.df$silhouette_category) / nrow(fa.sil.df)
# Third Step: Membership analysis: Values in text
# Membership analysis
famemb.fun <- function(i){
  df <- austria_cluster.df %>% filter(., id %in% i)
  fan = fanny(scale(df$Rating), k = 2)
  memshp <- fan$membership %>% as.data.frame(.)
  strong <- memshp %>% filter(., V1 < 0.25 | V1 > 0.75)
  strong.sh <- length(strong$V1)
  #  avg.mem <- mean(memshp$V1)
  return(strong.sh)
}

famemb.list <- lapply(unique(working.ids$id), FUN=famemb.fun)
unlist(famemb.list) / 6
famemb.graz <- data.frame(id = working.ids$id,
                      mmsh.value = unlist(famemb.list),
                      mmsh.share = unlist(famemb.list) / 6) %>%
  filter(., id %in% optclust2.df$id)
table(famemb.graz$mmsh.share)
####  END FANNY ################################################
################################################################





















## to be removed later:
### Robustness check: use fanny (fuzzy clustering) instead of kmeans
fanny_optimal_clusters <- function(df, id_var, rating_var, max_k = 5) {
  
  results <- df %>%
    group_by(!!sym(id_var)) %>%
    group_split() %>%
    lapply(function(sub_df) {
      id_val <- unique(sub_df[[id_var]])
      ratings <- sub_df[[rating_var]]
      # Remove duplicates and check how many unique values exist
      unique_ratings <- unique(ratings)
      num_unique <- length(unique_ratings)
      # Skip clustering if there are not enough unique values
      if (num_unique < 2) {
        return(data.frame(id = id_val, optimal_k = NA, best_silhouette = NA))
      }
      # Limit max_k to number of unique values
      max_k_adj <- 2 #min(max_k, num_unique)
      # Initialize silhouette scores
      silhouette_scores <- rep(NA, max_k_adj)
      
      for (k in 2:max_k_adj) {
        kmeans_result <- fanny(ratings, k = 2, metric ="euclidean")
        cluster_labels <- kmeans_result$cluster
        
        # Compute silhouette score only if at least two clusters exist
        if (length(unique(cluster_labels)) < 2) next
        
        sil_values <- silhouette(cluster_labels, dist(ratings))
        silhouette_scores[k] <- mean(sil_values[, 3], na.rm = TRUE)
      }
      
      # Determine the optimal number of clusters
      if (all(is.na(silhouette_scores))) {
        return(data.frame(id = id_val, optimal_k = NA, best_silhouette = NA))
      }
      
      optimal_k <- which.max(silhouette_scores[-1]) + 1
      best_silhouette <- max(silhouette_scores, na.rm = TRUE)
      
      return(data.frame(id = id_val, optimal_k = optimal_k, best_silhouette = best_silhouette))
    }) %>%
    bind_rows()
  
  return(results)
}
# Run the function on austria_cluster.df
fanny_clusters_df <- fanny_optimal_clusters(austria_cluster.df, "id", "Rating")

optclust.fanny <- fanny_clusters_df  %>%
  filter(id %in% optclust2.df$id) %>%
  mutate(
    silhouette_category = case_when(
      best_silhouette > 0.70 ~ "strong",
      best_silhouette > 0.50 ~ "moderate",
      best_silhouette > 0.25 ~ "weak",
      best_silhouette > 0    ~ "poor",
      best_silhouette <= 0   ~ "incorrect",
      is.na(best_silhouette) ~ NA_character_  # Keep NA if silhouette is missing
    )
  )
table(optclust.fanny$silhouette_category)
table(optclust.fanny$silhouette_category) / nrow(optclust2.df)

##


