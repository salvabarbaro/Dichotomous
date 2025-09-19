library(dplyr)
library(IC2)
library(dineq)

library(dplyr)
# library(IC2)  # ensure loaded

GiniComparison.fun <- function(i, data) {
  df <- data %>%
    filter(id == i) %>%
    filter(!is.na(Rating), !is.na(Approval)) %>%
    mutate(Approval = factor(Approval))

  # need both groups present
  if (nrow(df) == 0 || dplyr::n_distinct(df$Approval) < 2) return(NULL)

  # helper for safe extraction
  safe_extract <- function(obj, slot) {
    if (!is.null(obj$decomp[[slot]])) as.numeric(obj$decomp[[slot]]) else NA_real_
  }

  # --- YL decomposition
  yl <- try(decompSGini(x = df$Rating, z = df$Approval, decomp = "YL", ELMO = FALSE),
            silent = TRUE)

  yl.total   <- if (!inherits(yl, "try-error")) as.numeric(yl$ineq$index)     else NA_real_
  yl.within  <- if (!inherits(yl, "try-error")) as.numeric(yl$decomp$within)  else NA_real_
  yl.between <- if (!inherits(yl, "try-error")) as.numeric(yl$decomp$between) else NA_real_
  yl.overlap <- if (!inherits(yl, "try-error")) safe_extract(yl, "stratif")   else NA_real_

  # group-specific stratification contributions
  yl.stratif.0 <- if (!inherits(yl, "try-error") && !is.null(yl$stratif$contribStratifGroups["0"])) {
    as.numeric(yl$stratif$contribStratifGroups["0"])
  } else NA_real_

  yl.stratif.1 <- if (!inherits(yl, "try-error") && !is.null(yl$stratif$contribStratifGroups["1"])) {
    as.numeric(yl$stratif$contribStratifGroups["1"])
  } else NA_real_

  # --- BM decomposition
  bm <- try(decompSGini(x = df$Rating, z = df$Approval, decomp = "BM", ELMO = FALSE),
            silent = TRUE)

  bm.total   <- if (!inherits(bm, "try-error")) as.numeric(bm$ineq$index)     else NA_real_
  bm.within  <- if (!inherits(bm, "try-error")) as.numeric(bm$decomp$within)  else NA_real_
  bm.between <- if (!inherits(bm, "try-error")) as.numeric(bm$decomp$between) else NA_real_
  bm.overlap <- if (!inherits(bm, "try-error")) {
    if (!is.null(bm$decomp$overlap)) as.numeric(bm$decomp$overlap)
    else if (!is.null(bm$decomp$stratif)) as.numeric(bm$decomp$stratif)
    else NA_real_
  } else NA_real_

  # return tidy one-row df
  data.frame(
    id           = i,
    YL.tot       = yl.total,
    YL.within    = yl.within,
    YL.between   = yl.between,
    YL.overlap   = yl.overlap,
    YL.stratif.0 = yl.stratif.0,
    YL.stratif.1 = yl.stratif.1,
    BM.tot       = bm.total,
    BM.within    = bm.within,
    BM.between   = bm.between,
    BM.overlap   = bm.overlap
  )
}

# run and bind
ids <- unique(france_theil.df$id)
gcomp_list <- lapply(ids, GiniComparison.fun, data = france_theil.df)
gcomp <- dplyr::bind_rows(gcomp_list)  %>%
  mutate(WDP.YL = ifelse(YL.within < YL.between, 1, 0),
         WDP.BM = ifelse(BM.within < BM.between, 1, 0)        
        )
## Graz:
ids <- unique(austria_theil.df$id)
gcomp_list <- lapply(ids, GiniComparison.fun, data = austria_theil.df)
gcomp <- dplyr::bind_rows(gcomp_list)  %>%
  mutate(WDP.YL = ifelse(YL.within < YL.between, 1, 0),
         WDP.BM = ifelse(BM.within < BM.between, 1, 0)        
        )
# Grenoble
ids <- unique(grenoble_theil.df$id)
gcomp_list <- lapply(ids, GiniComparison.fun, data = grenoble_theil.df)
gcomp <- dplyr::bind_rows(gcomp_list)  %>%
  mutate(WDP.YL = ifelse(YL.within < YL.between, 1, 0),
         WDP.BM = ifelse(BM.within < BM.between, 1, 0)        
        )


#res0.yl <- gcomp %>% filter(., abs(YL.overlap) < 1e-3)  
#res0.bm <- gcomp %>% filter(., abs(BM.overlap) < 1e-3)  




ggplot(data = gcomp, aes(x = BM.within, y = YL.within)) +
  geom_point()

overlaps_long <- gcomp %>%
  select(id, BM.overlap, YL.overlap) %>%
  pivot_longer(cols = c(BM.overlap, YL.overlap),
               names_to = "method",
               values_to = "overlap")

# density plot
ggplot(overlaps_long, aes(x = overlap, fill = method, colour = method)) +
  geom_density(alpha = 0.3) +
  labs(x = "Overlap component", y = "Density",
       title = "Distribution of BM vs. YL overlap contributions") +
  theme_minimal(base_size = 22)

ggplot(overlaps_long, aes(x = method, y = overlap, fill = method)) +
  geom_boxplot(alpha = 0.6) +
  labs(x = "Method", y = "Overlap component",
       title = "BM vs. YL overlap contributions") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")

## within
within_long <- gcomp %>%
  select(id, BM.within, YL.within) %>%
  pivot_longer(cols = c(BM.within, YL.within),
               names_to = "method",
               values_to = "within")

ggplot(within_long, aes(x = within, fill = method, colour = method)) +
  geom_density(alpha = 0.3) +
  labs(x = "within component", y = "Density",
       title = "Distribution of BM vs. YL within contributions") +
  theme_minimal(base_size = 22)
#between
between_long <- gcomp %>%
  select(id, BM.between, YL.between) %>%
  pivot_longer(cols = c(BM.between, YL.between),
               names_to = "method",
               values_to = "between")

ggplot(between_long, aes(x = between, fill = method, colour = method)) +
  geom_density(alpha = 0.3) +
  labs(x = "between component", y = "Density",
       title = "Distribution of BM vs. YL between contributions") +
  theme_minimal(base_size = 22)



head(france_theil.df)
id0 <- france_theil.df %>% filter(., id == 4) %>% mutate(Approval = as.factor(Approval))
yl.0 <- decompSGini(x = id0$Rating, z = id0$Approval, decomp = "YL")
#yl.0$stratif
yl.0$decomp$stratif



########################################################################################
## path-independent within component!
library(ineq)
## use the code from Code_for_Barbaro.R
id0 <- france_theil.df %>% filter(., id == 0)
gini_within_add_path(incomes = id0$Rating, group = id0$Approval)


within.path.fun <- function(i, data){
  df <- data %>%
    filter(id == i) 
    if (nrow(df) == 0 || dplyr::n_distinct(df$Approval) < 2) return(NULL)
  bari <- gini_within_add_path(incomes = df$Rating, group = df$Approval)
#  return(bari[2])
  data.frame(id = i, giniwpath = bari[2])
}

ids <- unique(france_theil.df$id)
bari_list <- lapply(ids, within.path.fun, data = france_theil.df)
bari.df <- bari_list %>% bind_rows(.)

gcompBari <- gcomp %>% left_join(x = ., y = bari.df, by = "id") %>%
  dplyr::select(., c("id", "YL.within", "BM.within", "giniwpath"))

bari_long <- gcompBari %>%
  pivot_longer(cols = c(YL.within, BM.within, giniwpath),
               names_to = "method",
               values_to = "within")
ggplot(bari_long, aes(x = within, fill = method, colour = method)) +
  geom_density(alpha = 0.3) +
  labs(x = "within component", y = "Density",
       title = "Distribution of BM vs. YL vs. Path within contributions") +
  theme_minimal(base_size = 22)
