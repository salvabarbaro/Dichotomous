######################################################################################################
library(dplyr)
library(cluster)
library(factoextra)
library(ggplot2)
library(tidyr)
library(parallel)
library(modelsummary)

france22.df <- read.csv("DATA/France22.csv", header = T)
###
## Data summary (Supplementary Material)
modelsummary::datasummary(
  All(france22.df) ~ N + Mean + SD + Min + Median + Max,
  data = france22.df#,
#  output = "france22_summary.tex"
)

scale_to_range <- function(x, new_min, new_max) {
  old_min <- min(x, na.rm = T)
  old_max <- max(x, na.rm = T)
  scaled_x <- ((x - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
  return(scaled_x)
}

france_long.df <- france22.df %>%
  select(id, starts_with("AV_"), starts_with("EV_"))  %>%
  pivot_longer(cols = starts_with("EV_"), names_to = "Candidate", values_to = "Approval") %>%
  pivot_longer(cols = starts_with("AV_"), names_to = "Approval_Candidate", values_to = "Rating") %>%
  filter(substring(Candidate, 3) == substring(Approval_Candidate, 3)) %>%  # Ensure candidate names match 
  mutate(Rating = ifelse(Rating ==50, NA, Rating))    ## robustness check: convert 50 to NA

# Transformation... - apply function
newratings <- scale_to_range(x = france_long.df$Rating, 2,3)

france_theil.df <- france_long.df %>%
  mutate(Rating = newratings) 
rm(newratings)  ## clean up

france_cluster.df <- france_theil.df %>% 
  group_by(id) %>%
  filter(all(!is.na(Rating))) %>%
  ungroup()
###
rm(france_theil.df, france22.df, france_long.df, scale_to_range)
gc()
## we need: france_cluster.df
#
rowwise.kmeans.fun <- function(r, df){
  ro <- df %>% filter(id == r)
  res.clu <- kmeans(
    ro$Rating, centers = 2, nstart = 25
  )

  res.fit <- fitted(res.clu, method = "centers" )

  res2.fit <- fitted(res.clu, method = "classes")

  res.df <- data.frame(
    id = ro$id,
    Candidate = ro$Candidate,
    Approval = ro$Approval,
#    Approval_Candidate = ro$Approval_Candidate,
    Rating = ro$Rating,
    clus = unlist(res2.fit),
    value = unlist(res.fit)
  )
  return(res.df)
}

# only individuals with differences in their ratings
working.ids <- france_cluster.df %>%
  group_by(id) %>%
  reframe(l = var(Rating, na.rm = TRUE)) %>%
  filter(l > 0)

set.seed(55234)

france01 <- lapply(
  working.ids$id,
  rowwise.kmeans.fun,
  df = france_cluster.df
)

## start: cluster number transformation
## "1": approved, "2" disapproved. We transform this to 0: disapproved, 1: approved 
approv.cluster.fun <- function(df){

  df %>%
    mutate(
      clus.assign = ifelse(value < mean(value), 0, 1))
}

approv.cluster.df <- lapply(
  france01,
  approv.cluster.fun
) %>%
  bind_rows() %>%
  rename(Cluster = clus.assign) %>%
  mutate(match = ifelse(Approval == Cluster, 1, 0))

rownames(approv.cluster.df) <- seq_len(nrow(approv.cluster.df))

cat(
  "\nShare of approval/cluster-derived approval coincidences:\n"
)

print(
  table(approv.cluster.df$match) / nrow(approv.cluster.df)
)

rm(working.ids, france01, approv.cluster.fun, rowwise.kmeans.fun)
gc()


mod.match <- glm(
  match ~ Candidate,
  family = binomial(link = "logit"),
  data = approv.cluster.df
)

preds <- marginaleffects::avg_predictions(
  mod.match,
  by = "Candidate",
  vcov = ~id
)

preds


preds <- preds %>%
  mutate(
    Candidate = recode(
      Candidate,
      "EV_AH"  = "Hidalgo",  #10      centre
      "EV_EM"  = "Macron",  #1        centre
      "EV_EZ"  = "Zemmour", #4        extr right
      "EV_FR"  = "Roussel", #8        extr left
      "EV_JJ"  = "Lassalle",  #7?      NA 
      "EV_JLM" = "Mélenchon", #3      extr left
      "EV_MLP" = "Le Pen",  #2        extr right
      "EV_NA"  = "Arthaud", #12       extr left
      "EV_NDA" = "Dupont-Aignan", #9  extr right / centre
      "EV_PP"  = "Poutou",  #11       extr left
      "EV_VP"  = "Pecresse", #5       centre
      "EV_YJ"  = "Jadot"  #6          centre
    )
  )

##https://de.wikipedia.org/wiki/Pr%C3%A4sidentschaftswahl_in_Frankreich_2022

ggplot(preds,
       aes(x = reorder(Candidate, estimate),
           y = estimate,
           ymin = conf.low,
           ymax = conf.high)) +
  geom_pointrange(linewidth = 0.5) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Predicted probability of agreement"
  ) +
  theme_bw(base_size = 24)
#ggsave("matchfigFrance.pdf", width = 16, height = 9)

rm(preds, mod.match)
gc()

df <- approv.cluster.df %>%
  mutate(Approval = as.numeric(Approval)) %>%
  group_by(id) %>%
  reframe(
    sum.app = sum(Approval), 
    sum.clu = sum(Cluster)) %>%
  distinct()

###################################################################
## Core Question 6
## first: we need a list of opt_k for each id
nb.cores <- 18
k_max <- 6
nb_iter <- 25

set.seed(55234)

kmeans.id.france <- function(df_id, kmax = k_max, nstart = nb_iter) {
  df_id <- df_id %>%
    dplyr::filter(
      !is.na(Rating),
      is.finite(Rating)
    )

  ratings <- df_id$Rating

  # Mindestvoraussetzungen
  if (
    length(ratings) < 6 ||
    length(unique(ratings)) < kmax ||
    !is.finite(sd(ratings)) ||
    sd(ratings) == 0
  ) {
    return(NULL)
  }

  x <- matrix(ratings, ncol = 1)

  ks <- 2:kmax

  # Distance matrix
  d <- dist(x)

  fits <- lapply(ks, function(k) {

    km <- stats::kmeans(x, centers = k, nstart = nstart)

    sil <- mean(cluster::silhouette(km$cluster, d)[, "sil_width"])

    list(k = k, fit = km, silhouette = sil ) })

  sil_values <- vapply(fits, function(z) z$silhouette, numeric(1) )

  best <- which.max(sil_values)

  opt_k <- fits[[best]]$k
  final_fit <- fits[[best]]$fit

  df_id %>%
    mutate(
      opt_k = opt_k,
      cluster = final_fit$cluster,
      cluster_center =
        final_fit$centers[final_fit$cluster, 1]
    )
}

split.france <- split(
  france_cluster.df,
  france_cluster.df$id
)

future::plan(
  future::multisession,
  workers = nb.cores
)

res.france <- future.apply::future_lapply(
  split.france,
  kmeans.id.france,
  kmax = 6,
  nstart = 25,
  future.seed = TRUE
)

cluster.france.df <- dplyr::bind_rows(res.france) 

AE.df <- cluster.france.df %>%
  left_join(x = ., y = approv.cluster.df, by = c("id", "Candidate", "Approval", "Rating"))

approv.sums <- AE.df %>%
  group_by(id) %>%
  reframe(
    sum.app = sum(Approval), 
    sum.clu = sum(Cluster) ) %>%
  distinct()

AE.df <- AE.df %>% left_join(x = ., y = approv.sums, by = "id")

long <- AE.df %>%
  pivot_longer(
    c(sum.app, sum.clu),
    names_to = "type",
    values_to = "n_candidates"
  ) %>%
  mutate(
    type = factor(
      type,
      levels = c("sum.app", "sum.clu")
    ),
    optk = factor(opt_k, levels = 2:6)
  )

m <- lme4::lmer(
  n_candidates ~ type + factor(optk) + (1 | id),
  data = long
)
summary(m)

#### compatible ranking ballots
crb <- AE.df %>%
  group_by(id) %>%
  filter(
    any(Approval == 0 & !is.na(Rating)),
    any(Approval == 1 & !is.na(Rating))
  ) %>%
  summarise(
    maxDisappr = max(Rating[Approval == 0], na.rm = TRUE),
    minAppr    = min(Rating[Approval == 1], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    compatible.strong = ifelse(maxDisappr <  minAppr, 1, 0),
    compatible.weak   = ifelse(maxDisappr <= minAppr, 1, 0)
  )

head(crb)
sapply(crb, mean)
