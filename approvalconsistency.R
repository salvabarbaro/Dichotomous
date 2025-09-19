ggplot(data = grenoble_theil.df,
    aes(y = id, x = Rating,  group = factor(Approval), colour = factor(Approval))
  ) +
  geom_point()


threshold_consistency <- function(df) {
  
  # keep only ids with at least one approved and one disapproved candidate
  df_filtered <- df %>%
    group_by(id) %>%
    filter(any(Approval == 0 & !is.na(Rating)) &
           any(Approval == 1 & !is.na(Rating))) %>%
    ungroup()
  
  # total number of ids considered
  n_total <- n_distinct(df_filtered$id)
  
  # strict condition
  res_strict <- df_filtered %>%
    group_by(id) %>%
    summarise(cond = max(Rating[Approval == 0], na.rm = TRUE) >
                      min(Rating[Approval == 1], na.rm = TRUE),
              .groups = "drop")
  n_strict <- sum(res_strict$cond)
  
  # weak condition
  res_weak <- df_filtered %>%
    group_by(id) %>%
    summarise(cond = max(Rating[Approval == 0], na.rm = TRUE) >=
                      min(Rating[Approval == 1], na.rm = TRUE),
              .groups = "drop")
  n_weak <- sum(res_weak$cond)
  
  tibble(
    n_total       = n_total,
    n_strict      = n_strict,
    share_strict  = 1 - (n_strict / n_total),
    n_weak        = n_weak,
    share_weak    = 1 - (n_weak / n_total)
  )
}

threshold_consistency(df = france_theil.df)

## which id's do have strict threshold consistent preferences
strict_ids <- function(df) {
  df %>%
    group_by(id) %>%
    filter(any(Approval == 0 & !is.na(Rating)) &
             any(Approval == 1 & !is.na(Rating))) %>%
    summarise(cond = max(Rating[Approval == 0], na.rm = TRUE) >
                      min(Rating[Approval == 1], na.rm = TRUE),
              .groups = "drop") %>%
    filter(cond) %>%
    pull(id)
}

ids_strict <- strict_ids(france_theil.df)

###########################################################
## bootstrap confidence intervals
library(boot)
# Slightly modified to return only share_strict and share_weak
check_approval_consistency <- function(df) {
  
  df_filtered <- df %>%
    group_by(id) %>%
    filter(any(Approval == 0 & !is.na(Rating)) &
             any(Approval == 1 & !is.na(Rating))) %>%
    ungroup()
  
  n_total <- n_distinct(df_filtered$id)
  
  res_strict <- df_filtered %>%
    group_by(id) %>%
    summarise(cond = max(Rating[Approval == 0], na.rm = TRUE) >
                      min(Rating[Approval == 1], na.rm = TRUE),
              .groups = "drop")
  n_strict <- sum(res_strict$cond)
  
  res_weak <- df_filtered %>%
    group_by(id) %>%
    summarise(cond = max(Rating[Approval == 0], na.rm = TRUE) >=
                      min(Rating[Approval == 1], na.rm = TRUE),
              .groups = "drop")
  n_weak <- sum(res_weak$cond)
  
  return(c(
    share_strict = 1 - (n_strict / n_total),
    share_weak   = 1 - (n_weak / n_total)
  ))
}

# wrapper for boot at id level
boot_fun <- function(data, indices) {
  # data here will just be the vector of unique IDs
  ids <- data[indices]
  
  # reconstruct bootstrap sample
  d <- france_theil.df %>% filter(id %in% ids)  ## adjust data set!
  
  check_approval_consistency(d)
}

# number of unique respondents
ids <- unique(france_theil.df$id)   # adjust data set

set.seed(55234)
b_res <- boot(
  data = ids,            # <- note: pass vector of ids, not full df
  statistic = boot_fun,
  R = 1000,
  parallel = "multicore", # or "snow" on Windows
  ncpus = 8
)

# inspect
b_res$t0
boot.ci(b_res, type = "perc", index = 1) # CI for share_strict
boot.ci(b_res, type = "perc", index = 2) # CI for share_weak

###################################################################
#### Midpoint calculation
calc_midpoints <- function(df) {
  df %>%
    group_by(id) %>%
    filter(any(Approval == 0 & !is.na(Rating)) &
             any(Approval == 1 & !is.na(Rating))) %>%
    summarise(
      min_approved = min(Rating[Approval == 1], na.rm = TRUE),
      max_disapproved = max(Rating[Approval == 0], na.rm = TRUE),
      midpoint = (min_approved + max_disapproved) / 2,
      median.rating = median(Rating, na.rm = T),
      .groups = "drop"
    )
}

midpGrenoble   <- calc_midpoints(grenoble_theil.df)
midpGraz       <- calc_midpoints(austria_theil.df)
midpFrance22NA <- calc_midpoints(france_theil.df)
midpFrance22  <-  calc_midpoints(france_theil.df)

m1 <- midpGrenoble   %>% mutate(dataset ="Grenoble")
m2 <- midpGraz       %>% mutate(dataset = "Graz")
m3 <- midpFrance22   %>% mutate(dataset = "France22")
m4 <- midpFrance22NA %>% mutate(dataset = "France22NA")
all_mid <- bind_rows(m1, m2, m3, m4)

all_mid <- all_mid %>%
  mutate(dataset = factor(dataset,
                          levels = c("Grenoble", "Graz", "France22", "France22NA")))

ggplot(all_mid, aes(x = dataset, y = midpoint)) +
  geom_violin(trim = FALSE, alpha = 0.6, fill = "steelblue", col ="darkred") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, col = "darkred") +
  labs(y = "Threshold", x = "") + 
  theme(legend.position = "none") +
  theme_gray(base_size = 22)
ggsave("~/Documents/Research/Dichotomous/git/67b5f34c104b85acf4a11317/thresholdDist.pdf", width = 16, height = 9)


m1 %>% mutate(mAm = ifelse(midpoint > median.rating, 1, 0)) %>% summarize(l = length(mAm), l.sh = sum(mAm)/nrow(.)) ## 91.4%
m2 %>% mutate(mAm = ifelse(midpoint > median.rating, 1, 0)) %>% summarize(l = length(mAm), l.sh = sum(mAm)/nrow(.)) ## 65.6%
m3 %>% mutate(mAm = ifelse(midpoint > median.rating, 1, 0)) %>% summarize(l = length(mAm), l.sh = sum(mAm)/nrow(.)) ## 88.7%
m4 %>% mutate(mAm = ifelse(midpoint > median.rating, 1, 0)) %>% summarize(l = length(mAm), l.sh = sum(mAm)/nrow(.)) ## 90.9%
