nb.cores <- 18

## Define the max number of clusters
k_max <- 4

pam.id <- function(df_id, kmax = k_max) {

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

  # distance matrix
  d <- dist(x)

  fits <- lapply(ks, function(k) {

    pm <- cluster::pam(
      x,
      k = k,
      metric = "euclidean"
    )

    sil <- mean(
      cluster::silhouette(
        pm$clustering,
        d
      )[, "sil_width"]
    )

    list(
      k = k,
      fit = pm,
      silhouette = sil
    )
  })

  sil_values <- vapply(
    fits,
    function(z) z$silhouette,
    numeric(1)
  )

  best <- which.max(sil_values)

  opt_k <- fits[[best]]$k
  final_fit <- fits[[best]]$fit

  # PAM-Medoids
  medoids <- as.numeric(final_fit$medoids)

  df_id %>%
    mutate(
      opt_k = opt_k,
      cluster = final_fit$clustering,
      cluster_medoid =
        medoids[final_fit$clustering]
    )
}

split.df <- split(
  cs_p.df,
  interaction(
    cs_p.df$case_ID,
    cs_p.df$id,
    drop = TRUE
  )
)

res.pam <- parallel::mclapply(
  split.df,
  pam.id,
  kmax = 4,
  mc.cores = nb.cores
)

cluster.pam.df <- dplyr::bind_rows(res.pam)

kmeans.short <- cluster.df %>%
  dplyr::select(case_ID, id, opt_k) %>%
  distinct() %>%
  rename(opt_k_kmeans = opt_k)

pam.short <- cluster.pam.df %>%
  dplyr::select(case_ID, id, opt_k) %>%
  distinct() %>%
  rename(opt_k_pam = opt_k)

compare.k <- left_join(
  kmeans.short,
  pam.short,
  by = c("case_ID", "id")
)

pam.short <- cluster.pam.df %>%
  select(case_ID, id, opt_k) %>%
  distinct()

pam.table <- pam.short %>%
  mutate(
    opt_k_cat = ifelse(opt_k >= 4, 4, opt_k)
  ) %>%
  count(case_ID, opt_k_cat) %>%
  pivot_wider(
    names_from = opt_k_cat,
    values_from = n,
    values_fill = 0,
    names_prefix = "k"
  ) %>%
  mutate(
    total = k2 + k3 + k4,
    k2_pct = 100 * k2 / total,
    k3_pct = 100 * k3 / total,
    k4_pct = 100 * k4 / total
  )


common.df <- left_join(
  cou.table %>%
    select(case_ID, k2_pct) %>%
    rename(kmeans = k2_pct),
  pam.table %>%
    select(case_ID, k2_pct) %>%
    rename(pam = k2_pct),
  by = "case_ID"
)

cor.test(common.df$kmeans,
         common.df$pam)

#Pearson's product-moment correlation
#
#data:  common.df$kmeans and common.df$pam
#t = 105.32, df = 171, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9897206 0.9943527
#sample estimates:
#      cor 
#0.9923796 
ggplot(common.df, aes(x = kmeans, y = pam)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
#  geom_jitter(width = .05, height = .05) +
  geom_point(size = 2.5, alpha = .8) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "steelblue"
  ) +
  coord_equal() +
  labs(
    x = "Share with optimal k = 2 (k-means)",
    y = "Share with optimal k = 2 (PAM)"
  ) +
  theme_bw(base_size = 20)
ggsave("PAMkmeansCompare.pdf", width = 16, height = 9)


##### test