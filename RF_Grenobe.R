setwd("~/Documents/Research/Dichotomous")

library(cluster)
library(dplyr)
library(ggplot2)
library(factoextra)
library(boot)
library(parallel)
library(tidyr)
library(IC2)

######################################################################################################

## Data preparation

grenoble.df <- read.csv("Data/GrenobleData_2017_Presid+.csv", sep = ";", header = T) %>%
  filter(!if_all(starts_with("EV_"), is.na)) %>%
  mutate(across(starts_with("EV_"), ~ na_if(., " None"))) %>%
  mutate(across(starts_with("EV_"), ~ as.numeric(.)) ) %>%
  mutate(across(starts_with("EV_"), ~ ifelse(. < 0 | . > 1, NA, .))) %>%
  mutate(across(c("AGE", "GENDER", "EDUC", "WORK"), ~ na_if(., "NSPP"))) %>%
  rename(id = VOTER)

grenoble_theil.df <- grenoble.df %>%
  rename(temp1 = AV_OPINION, temp2 = EV_OPINION) %>%
  select(id, starts_with("AV_"), starts_with("EV_"))  %>%
  pivot_longer(cols = starts_with("EV_"), names_to = "Candidate", values_to = "Rating") %>%
  pivot_longer(cols = starts_with("AV_"), names_to = "Approval_Candidate", values_to = "Approval") %>%
  mutate(Rating = Rating + 2) %>%
  filter(substring(Candidate, 3) == substring(Approval_Candidate, 3))  # Ensure candidate names match

grenoble_cluster.df <- grenoble_theil.df %>% 
  group_by(id) %>%
  filter(all(!is.na(Rating))) %>%
  ungroup()

#######################################################
### Section 3: Decomposition analysis
ic2decomp.fun <- function(i) {
  tryCatch({
    df <- grenoble_theil.df %>% 
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

ic2res <- lapply(unique(grenoble_theil.df$id), ic2decomp.fun)
ic2res.df <- ic2res %>% do.call(rbind, .) %>%
  mutate(WDP.Theil = ifelse(Theil.within < Theil.between, 1, 0),
         WDP.Gini =  ifelse(Gini.within < Gini.between, 1, 0),
         WDP.Atkinson = ifelse(Atkinson.within < Atkinson.between, 1, 0))
head(ic2res.df)

## Function to calculate the respective share of respondents with WDP
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
                                           
##########################################################################
### SECTION 4: Cluster Analysis
## Cluster analysis
# Step 1: We remove all id's without variance in the Rating
# Function working.ids selects the appropriate id-values 
working.ids <- grenoble_cluster.df %>% 
  group_by(id) %>% 
  reframe(l = var(Rating, na.rm = T)) %>%
  filter(., l > 0) 
kmax.value <- 6  # 5for Graz, 10 for Grenoble

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
                          df = grenoble_cluster.df, mc.cores = 14)

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
sil.values <- lapply(k2.ids$id, sil.fun, df = grenoble_cluster.df)
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
kmax.value = 4 # 
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
                             df = grenoble_cluster.df, mc.cores = 14)
optclust_fa.df <- data.frame(id = working.ids$id, 
                             optk = unlist(optclust_fa.list))
table(optclust_fa.df$optk[is.na(optclust_fa.df$optk)==F])
# The following values in Table
table(optclust_fa.df$optk[is.na(optclust_fa.df$optk)==F])/nrow(optclust_fa.df)

## Second step: Evaluate the silhouette width for each respondent. 
## Values in Table!
fasil.fun <- function(i){
  df <- grenoble_cluster.df %>% filter(., id %in% i)
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
  df <- grenoble_cluster.df %>% filter(., id %in% i)
  fan = fanny(scale(df$Rating), k = 2)
  memshp <- fan$membership %>% as.data.frame(.)
  strong <- memshp %>% filter(., V1 < 0.25 | V1 > 0.75)
  strong.sh <- length(strong$V1)
  #  avg.mem <- mean(memshp$V1)
  return(strong.sh)
}

famemb.list <- lapply(unique(working.ids$id), FUN=famemb.fun)
unlist(famemb.list) / 11
famemb.grenoble <- data.frame(id = working.ids$id,
                          mmsh.value = unlist(famemb.list),
                          mmsh.share = unlist(famemb.list) / 11) %>%
  filter(., id %in% optclust2.df$id)
table(famemb.grenoble$mmsh.share)
####  END FANNY ################################################
################################################################
