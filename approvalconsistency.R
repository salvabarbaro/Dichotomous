ggplot(data = grenoble_theil.df,
    aes(y = id, x = Rating,  group = factor(Approval), colour = factor(Approval))
  ) +
  geom_point()


check_approval_consistency <- function(df) {
  
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
    share_strict  = n_strict / n_total,
    n_weak        = n_weak,
    share_weak    = n_weak / n_total
  )
}

check_approval_consistency(df = grenoble_theil.df)
