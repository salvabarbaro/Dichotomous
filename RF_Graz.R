### Graz Dataset
### I. Read Data and prepare data
setwd("~/Documents/Research/Dichotomous/github/Dichotomous")

#########################################
library(dplyr)
library(cluster)
library(factoextra)
library(ggplot2)
library(ineq)
library(dineq)
library(rstatix)
library(tidyr)
library(IC2)  # use remotes::install_version("IC2"), the library is no longer maintained.
library(parallel)
library(gtsummary)
library(modelsummary)
library(rlang)
library(readxl)
#########################################
# Section 0: Read data and data handling
austria.df <- readxl::read_excel("DATA/Steirische_LTW_2019_Daten_Barbaro.xlsx") %>%
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

# extremistic
extr.df <- austria.df %>% 
  mutate(ext.right = ifelse(Rat.FPÖ > 17.5, 1, 0),
         ext.left =  ifelse(Rat.KPÖ > 17.5, 1, 0)) %>%
  select(., c("id", "ext.left", "ext.right"))
write.csv(extr.df, "extrGraz.csv", row.names = F)
rm(extr.df)
# Data frame for Theil analysis
austria_long.df <- austria.df %>%
  select(id, starts_with("Approv."), starts_with("Rat."))  %>%
  pivot_longer(cols = starts_with("Rat."), names_to = "Party", values_to = "Rating") %>%
  pivot_longer(cols = starts_with("Approv."), names_to = "Approval_Party", values_to = "Approval") %>%
  filter(substring(Party, 5) == substring(Approval_Party, 8)) %>%  # Ensure party names match
  mutate(Party = gsub("^Rat\\.", "", Party)) %>%
  select(., -c("Approval_Party"))

## extremist
#extrem.fun <- function(i){
#  t1 <- austria_long.df %>% filter(., id == i)
#  t2 <- t1$Party[which.max(t1$Rating)]
#  
#  return(t2)
#}
#plurality.df <- lapply(unique(austria_long.df$id), extrem.fun) %>% unlist(.)
#austria2.df <- austria.df %>% mutate(plurality = plurality.df)#

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
rm(newratings)
# Data frame for cluster analysis (na_omit of Rating)
austria_cluster.df <- austria_theil.df %>% 
  group_by(id) %>%
  filter(all(!is.na(Rating))) %>%
  ungroup()

#############################################################
### Section 3: Decomposition analysis
# Running the function using lapply, removing NULL results
load("AuxFunctions.RData")

ic2res <- lapply(unique(austria_theil.df$id), ic2decomp.fun, data = austria_theil.df)
ic2res.df <- ic2res %>% do.call(rbind, .) %>%
  mutate(WDP.Theil = ifelse(Theil.within < Theil.between, 1, 0),
         WDP.Gini =  ifelse(Gini.within < Gini.between, 1, 0),
         WDP.Atkinson = ifelse(Atkinson.within < Atkinson.between, 1, 0),
         WDP.SCV = ifelse(SCV.within < SCV.between, 1, 0)
        )
head(ic2res.df)

## Store ID's with WDP in reg.austria [for Secion on Sociodemographics]
reg.austria <- austria.df %>% 
  select(., c("id", "Gender", "Age", "Educ" )) %>%
  left_join(x = ., 
            y = ic2res.df %>% select(., c("id", "WDP.Theil", "WDP.Gini", "WDP.Atkinson")),
            by = "id")
# write.csv(reg.austria, "DATA/regaustria.csv", row.names = F)
##################################################################

## Values for Table 2:
compute_wdp_shares(ic2res.df)

####################################################################
## Bootstrap
cluster_bootstrap <- function(df, id_col, R = 1000) {
  unique_ids <- unique(df[[id_col]])  # Unique id values
  boot_results <- matrix(NA, nrow = R, ncol = 4)  # Store bootstrap results
  
  for (r in 1:R) {
    sampled_ids <- sample(unique_ids, replace = TRUE)  # Resample IDs
    boot_df <- df %>% filter(id %in% sampled_ids)  # Extract corresponding rows
    
    boot_results[r, ] <- compute_wdp_shares(boot_df)  # Compute WDP shares
  }
  
  colnames(boot_results) <- c("Theil", "Gini", "Atkinson", "SCV")
  return(as.data.frame(boot_results))
}

# Running the clustered bootstrap
set.seed(55234)  
boot_results <- cluster_bootstrap(ic2res.df, id_col = "id", R = 1000)
# Compute confidence intervals by percentile method
apply(boot_results, 2, quantile, probs = c(0.025, 0.975))  # 95% CI
##########################################################################
## Robustness Check according to the proposed method by Fleurbaey, Lambert,... 
out1 <- fleurbaey_pipeline(austria_theil.df)

# With your external table ic2res.df (must have columns id and Gini.between):
out2 <- fleurbaey_pipeline(
  data = austria_theil.df,
  join_df = ic2res.df,
  id_col = "id",
  income_col = "Rating",
  group_col = "Approval",
  join_id_col = "id",
  join_between_col = "BM.between"
)

#out2$res_ids
#out2$fleurbaey
out2$wdp_table
#out2$resnew_summ
rm(boot_results, ic2res, ic2res.df, bm_between_vec, cluster_bootstrap, compute_wdp_shares, cv_decomp_bm, 
   fleurbaey_pipeline, gini_between_add, gini_coefficient, gini_within_add_path)
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
strong.ids.graz <- silh.df$id[silh.df$silhouette_category=="strong"]
modera.ids.graz <- silh.df$id[silh.df$silhouette_category !="weak"]
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
#### SECTION 6: REASONABLE ASSIGNEMENT
## Step 1: we find out which candidate is optimally assigned to which cluster
##    The information is in kmeans()$classes
rowwise.kmeans.fun <- function(r, df){
  ro <- df %>% filter(., id %in% r) 
  res.clu <- kmeans(ro$Rating, centers = 2, nstart = 10)
  res.fit <- fitted(res.clu, method = "centers")  
  res2.fit <- fitted(res.clu, method = "classes")
  res.df <- data.frame(clus = unlist(res2.fit),
                       value = unlist(res.fit))
  # return(res.fit)
  return(res.df)
}

## select id's
working.ids <- austria_cluster.df %>% 
  group_by(id) %>% 
  reframe(l = var(Rating, na.rm = T)) %>%
  filter(., l > 0)

# by graz01, we obtain a list: the clustering (1 or 2) and the rating values
graz01 <- lapply(unique(working.ids$id), 
                 FUN = rowwise.kmeans.fun, 
                 df = austria_cluster.df)

## start: cluster number transformation
## "1": approved, "2" disapproved. We transform this to 0: disapproved, 1: approved 
approv.cluster.fun <- function(i){
  df <- graz01[i] %>% 
    as.data.frame(.) %>% 
    mutate(clus.assing = ifelse(value < mean(value), 0, 1))
}

approv.cluster.df <- lapply(X =1:nrow(working.ids), 
                            FUN = approv.cluster.fun) %>%
  do.call(rbind, .)
## end: cluster number transformation
##
## graz.df: a data frame indicating who was actually approved and who 'should' have been approved
graz.df <- cbind(austria_cluster.df %>% filter(., id %in% working.ids$id),
                 approv.cluster.df) %>%
  mutate(match = ifelse(Approval == clus.assing, 1, 0)) %>%
  left_join(x = .,
            y = austria.df %>% select(., c("id", "Gender", "Age", "Educ")),
            by = "id") 

## A. Similarities: how many parties are identical?
sum(graz.df$match[is.na(graz.df$match) == F]) / 
  length(graz.df$match[is.na(graz.df$match) == F])
# B. Match by party
parties <- unique(graz.df$Party)
match.fun <- function(party, df){
  df2 <- df %>% filter(., Party == party)
  match.sh <- sum(df2$match[is.na(df2$match) == F]) / 
    length(df2$match[is.na(df2$match) == F])
}
partymatch.list <- lapply(parties, match.fun, df = graz.df) 
partymatch.df <-   data.frame(
  Party = parties, 
  match.share = unlist(partymatch.list))

## Bootstrap to calculate CI
boot_match.fun <- function(data, indices, cand) {
  df_boot <- data[indices, ]  
  match.fun(cand, df_boot)    
}

bootstrap_results <- lapply(parties, function(cand) {
  boot.out <- boot(data = graz.df, 
                   statistic = function(d, i) boot_match.fun(d, i, cand), 
                   R = 1000, 
                   parallel = "multicore",
                   ncpus = 14)
  
  # Extract mean and confidence intervals (percentile)
  ci <- boot.ci(boot.out, type = "perc")$percent[4:5]
  
  # Return a named list
  list(
    Party = cand,
    Match_Share = mean(boot.out$t),
    Lower_CI = ci[1],
    Upper_CI = ci[2]
  )
})

candmatch.df <- do.call(rbind, lapply(bootstrap_results, as.data.frame))
candmatch.df
#stargazer::stargazer(candmatch.df, summary = F, rownames = F)
## Hypothesis: Individuals with \tilde{k}=2 exhibit a higher matching rate than 
## individuals with k!= 2. 
# 1. 
#phi coefficient of correlation between two dichotomous variables
phi.fun <- function(i, df) {
  tryCatch({
    data <- df %>% filter(id %in% i) %>% select(Approval, clus.assing)
    pp <- psych::phi(table(data))
    return(pp)
  }, error = function(e) {
    message("Error for id ", i, ": ", e$message)
    return(NA)
  })
}

phi.list <- lapply(unique(graz.df$id), phi.fun, df = graz.df)
#summary(unlist(phi.list) )
#phi.df <- data.frame(id = unique(greno.df$id), 
#                     phi.value = unlist(phi.list))
educ.df <- data.frame(Educ = unique(graz.df$Educ),
                      Educ.lvl = c(6, 5, 2, 3, 4, 1, NA))
age.df <- data.frame(Age = sort(unique(graz.df$Age)),
                     Age.num = 1:7)
graz.avg <- graz.df %>% group_by(id) %>%
  reframe(avg.match = mean(match, na.rm = T)) %>%
  mutate(optk2 = ifelse(id %in% optclust2.df$id, "k=2", "k=2+")) %>%
  mutate(phi.value = unlist(phi.list)) %>%
  left_join(x = ., 
            y = austria.df %>% 
              select(., c("id", "Age", "Gender", "Educ")),
            by = "id") %>% 
  left_join(x = ., y = educ.df, by = "Educ") %>%
  left_join(x = ., y = age.df, by = "Age")

#ols02 <- lm(formula = phi.value ~ optk2 + Age.num + Gender + Educ.lvl,
#            data = graz.avg)
#broom::tidy(logreg01, conf.int = T)
#summary(ols02)
#gtsummary::tbl_regression(ols02, conf.level = 0.95)
#texreg::texreg(ols02, single.row = T, label = "tb.ols", booktabs = T,
#               custom.model.names = "Graz")
write.csv(graz.avg, "Data/OLSgraz.csv", row.names = F)





## Consider socio-demographics
library(lme4)  # for random-intercept models  (to consider id)
educ.df <- data.frame(Educ = unique(graz.df$Educ),
                      Educ.lvl = c(6, 5, 2, 3, 4, 1, NA))
graz.df <- graz.df %>% left_join(x = .,
                                 y = educ.df, 
                                 by = "Educ")

logreg01 <- glmer(match ~ Party + Gender + 
                    as.numeric(Age) + Educ.lvl + 
                    (1 | id), 
                  family = binomial, data = graz.df)
#summary(logreg01)
