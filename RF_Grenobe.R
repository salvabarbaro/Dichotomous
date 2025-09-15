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
######################################################################################################
## Data preparation
grenoble.df <- read.csv("DATA/GrenobleData_2017_Presid+.csv", 
                        sep = ";", header = TRUE) %>%
  select(-c("AV_OPINION", "EV_OPINION")) %>%
  filter(!if_all(starts_with("EV_"), is.na)) %>%
  mutate(across(starts_with("EV_"), ~ na_if(., " None"))) %>%
  mutate(across(starts_with("EV_"), ~ gsub(",", ".", .) %>% trimws())) %>%  # Convert commas & trim spaces
  mutate(across(starts_with("EV_"), ~ ifelse(grepl("^10[0-6]\\.", .), NA, .))) %>%  # Set "100.", "101.", ..., "106." to NA
  mutate(across(starts_with("EV_"), ~ ifelse(grepl("^10\\.", .), NA, .))) %>%
  mutate(across(starts_with("EV_"), as.numeric)) %>%  # Convert to numeric
  mutate(across(starts_with("EV_"), ~ ifelse(. < 0 | . > 1, NA, .))) %>%
  mutate(across(c("AGE", "GENDER", "EDUC", "WORK"), ~ na_if(., " NSPP"))) %>%
  rename(id = VOTER)

ext.rights <- c(" NDA", " MLP")
ext.left   <- c(" NA" , " PP", " JLM")
ids.radright_gre <- grenoble.df$id[grenoble.df$OFFIC %in% ext.rights]
ids.radleft_gre  <- grenoble.df$id[grenoble.df$OFFIC %in% ext.left]
#write.csv(ids.radright_gre, "idsradrightgre.csv", row.names = F)
#write.csv(ids.radleft_gre,  "idsradleftgre.csv", row.names = F)
rm(ids.radleft_gre, ids.radright_gre, ext.left, ext.rights)

grenoble_theil.df <- grenoble.df %>%
#  rename(temp1 = AV_OPINION, temp2 = EV_OPINION) %>%
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
load("AuxFunctions.RData")

ic2res <- mclapply(unique(grenoble_theil.df$id), ic2decomp.fun, data = grenoble_theil.df, mc.cores = 8)
ic2res.df <- ic2res %>% do.call(rbind, .) %>%
  mutate(WDP.Theil = ifelse(Theil.within < Theil.between, 1, 0),
         WDP.Gini =  ifelse(Gini.within < Gini.between, 1, 0),
         WDP.Atkinson = ifelse(Atkinson.within < Atkinson.between, 1, 0),
         WDP.SCV = ifelse(SCV.within < SCV.between, 1, 0)
        )
head(ic2res.df)

## Store ID's with WDP in reg.grenoble [for sec 7]
reg.grenoble <- grenoble.df %>% 
  select(., c("id", "GENDER", "AGE", "EDUC")) %>%
  left_join(x = ., 
            y = ic2res.df %>% select(., c("id", "WDP.Theil", "WDP.Gini", "WDP.Atkinson")),
            by = "id")
#write.csv(reg.grenoble, "reggrenoble.csv", row.names = F)

## Values for Table 2:
compute_wdp_shares(ic2res.df)
####################################################################
## Bootstrap
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
out1 <- fleurbaey_pipeline(grenoble_theil.df)

# With your external table ic2res.df (must have columns id and Gini.between):
out2 <- fleurbaey_pipeline(
  data = grenoble_theil.df,
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
## for section 7:
strong.ids.grenoble <- silh.df$id[silh.df$silhouette_category=="strong"]
modera.ids.grenoble <- silh.df$id[silh.df$silhouette_category !="weak"]
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
working.ids <- grenoble_cluster.df %>% 
  group_by(id) %>% 
  reframe(l = var(Rating, na.rm = T)) %>%
  filter(., l > 0)

# by grenoble01, we obtain a list: the clustering (1 or 2) and the rating values
grenoble01 <- lapply(unique(working.ids$id), 
                 FUN = rowwise.kmeans.fun, 
                 df = grenoble_cluster.df)

## start: cluster number transformation
## "1": approved, "2" disapproved. We transform this to 0: disapproved, 1: approved 
approv.cluster.fun <- function(i){
  df <- grenoble01[i] %>% 
    as.data.frame(.) %>% 
    mutate(clus.assing = ifelse(value < mean(value), 0, 1))
}

approv.cluster.df <- lapply(X =1:nrow(working.ids), 
                            FUN = approv.cluster.fun) %>%
  do.call(rbind, .)
## end: cluster number transformation
##
## greno.df: a data frame indicating who was actually approved and who 'should' have been approved
greno.df <- cbind(grenoble_cluster.df %>% 
                       filter(., id %in% working.ids$id),
                 approv.cluster.df) %>%
  mutate(match = ifelse(Approval == clus.assing, 1, 0)) %>%
  left_join(x = .,
            y = grenoble.df %>% select(., c("id", "AGE", "GENDER", "EDUC", "WORK")),
            by = "id")  %>%
  mutate(GENDER = na_if(GENDER, " NSPP"),
         EDUC =   na_if(EDUC,   " NSPP")) %>%
  mutate(Educ.lvl = as.numeric(ifelse(EDUC  == " S", 3, EDUC)))

## A. Similarities: how many parties are identical?
sum(greno.df$match[is.na(greno.df$match) == F]) / 
  length(greno.df$match[is.na(greno.df$match) == F])
# B. Match by party
candidates <- unique(greno.df$Candidate)
match.fun <- function(cand, df){
  df2 <- df %>% filter(., Candidate == cand)
  match.sh <- sum(df2$match[is.na(df2$match) == F]) / 
    length(df2$match[is.na(df2$match) == F])
  return(match.sh)
}
candmatch.list <- lapply(candidates, match.fun, df = greno.df) 
candmatch.df <-   data.frame(
  Candidates = candidates, 
  match.share = unlist(candmatch.list))
######################################################
## Bootstrap to calculate CI
boot_match.fun <- function(data, indices, cand) {
  df_boot <- data[indices, ]  
  match.fun(cand, df_boot)    
}

bootstrap_results <- lapply(candidates, function(cand) {
  boot.out <- boot(data = greno.df, 
                   statistic = function(d, i) boot_match.fun(d, i, cand), 
                   R = 1000, 
                   parallel = "multicore",
                   ncpus = 14)
  
  # Extract mean and confidence intervals (percentile)
  ci <- boot.ci(boot.out, type = "perc")$percent[4:5]
  
  # Return a named list
  list(
    Candidate = cand,
    Match_Share = mean(boot.out$t),
    Lower_CI = ci[1],
    Upper_CI = ci[2]
  )
})

candmatch.df <- do.call(rbind, lapply(bootstrap_results, as.data.frame))
candmatch.df
stargazer::stargazer(candmatch.df, summary = F, rownames = F)
##############################################
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

phi.list <- lapply(unique(greno.df$id), phi.fun, df = greno.df)
#summary(unlist(phi.list) )
#phi.df <- data.frame(id = unique(greno.df$id), 
#                     phi.value = unlist(phi.list))

greno.avg <- greno.df %>% group_by(id) %>%
  reframe(avg.match = mean(match, na.rm = T)) %>%
  mutate(optk2 = ifelse(id %in% optclust2.df$id, "k=2", "k=2+")) %>%
  mutate(phi.value = unlist(phi.list)) %>%
  left_join(x = ., 
            y = grenoble.df %>% 
              select(., c("id", "AGE", "GENDER", "EDUC", "WORK")),
            by = "id") %>%
  mutate(Educ.lvl = as.numeric(ifelse(EDUC  == " S", 3, EDUC)))

ols01 <- lm(formula = phi.value ~ optk2 + as.numeric(AGE) + GENDER + Educ.lvl,
                data = greno.avg)
#broom::tidy(logreg01, conf.int = T)
summary(ols01)
gtsummary::tbl_regression(ols01)
texreg::texreg(ols01, single.row = T, label = "tb.ols", booktabs = T,
               custom.model.names = "Grenoble")
#write.csv(greno.avg, "Data/OLSgrenoble.csv", row.names = F)

### Result: persons with opt k > 2 exhibit a significant lower phi coefficient. 


#boot.pval::boot_summary(lm(formula = phi.value ~ optk2  + GENDER + Educ.lvl,
#                           data = greno.avg))
# ## Compare \tilde{k}=2 and >2
# greno2.df <- greno.df %>%
#   mutate(optk2 = ifelse(id %in% optclust2.df$id, "k=2", "k=2+"))
# #  adjustment
# match2.fun <- function(cand, df, g){
#   df2 <- df %>% filter(., Candidate == cand, optk2 == g)
#   match.sh <- sum(df2$match[is.na(df2$match) == F]) / 
#     length(df2$match[is.na(df2$match) == F])
#   return(match.sh)
# }
# #
# candmatch.listk2 <-   lapply(candidates, match2.fun, df = greno2.df, g = "k=2") 
# candmatchk2.df <-   data.frame(
#   Candidates = candidates, 
#   match.share = unlist(candmatch.listk2))
# ##
# candmatch.listk3 <-   lapply(candidates, match2.fun, df = greno2.df, g = "k=2+") 
# candmatchk3.df <-   data.frame(
#   Candidates = candidates, 
#   match.share = unlist(candmatch.listk3))
# candmatch.comp <- candmatchk2.df %>% 
#   left_join(x = ., 
#             y = candmatchk3.df, 
#             by = "Candidates") %>% 
#   setNames(c("Candidates", "Match.k2", "Match.k3"))

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
                  family = binomial, data = grenoble.df)
#summary(logreg01)
#broom.mixed::tidy(logreg01, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
#texreg::screenreg(logreg01)
gtsummary::tbl_regression(logreg01, exponentiate = T)
## same model without random effect (similar effects)
logreg02 <- glm(match ~ Party + Gender + 
                  as.numeric(Age) + Educ.lvl,
                family = "binomial", data = graz.df)
gtsummary::tbl_regression(logreg02, exponentiate = T)
## model without parties
logreg03 <- glm(match ~ Gender + 
                  as.numeric(Age) + Educ.lvl,
                family = "binomial", data = graz.df)
gtsummary::tbl_regression(logreg03, exponentiate = T)
# model with parties as random effect
logreg04 <- glmer(match ~ Gender + 
                    as.numeric(Age) + Educ.lvl + 
                    (1 | Party), 
                  family = binomial, data = graz.df)
gtsummary::tbl_regression(logreg04, exponentiate = T)
