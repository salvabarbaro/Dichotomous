# Gini helper
gini <- function(x) {
  n <- length(x)
  if (n <= 1) return(0)
  xs <- sort(x)
  (2 * sum(xs * seq_len(n)) / (n * sum(xs))) - (n + 1) / n
}

# Bhattacharya–Mahalanobis subgroup decomposition of the Gini
# GBM_w = sum_k (p_k^2 * (mu_k / mu) * G_k)
# GBM_b = (1 / (2*mu)) * sum_k sum_h p_k p_h * |mu_k - mu_h|
bm_gini_bm <- function(y, g) {
  N  <- length(y)
  mu <- mean(y)

  by_mean <- tapply(y, g, mean)
  by_size <- tapply(y, g, length)
  by_gini <- tapply(y, g, gini)
  p       <- by_size / N

  # Within (BM)
  within <- sum((p^2) * (by_mean / mu) * by_gini)

  # Between (BM)
  # double sum over groups of p_k p_h |mu_k - mu_h|, then divide by 2*mu
  M <- expand.grid(k = seq_along(by_mean), h = seq_along(by_mean))
  between <- sum(p[M$k] * p[M$h] * abs(by_mean[M$k] - by_mean[M$h])) / (2 * mu)

  overall <- gini(y)
  overlap <- overall - (within + between)  # 0 if no overlap

  list(
    overall_gini = overall,
    within_bm    = within,
    between_bm   = between,
    overlap      = overlap,
    check_sum    = within + between
  )
}

# ----- Your example -----
df <- tibble::tibble(
  Approval = c(0,0,0,0,1,0,1,0,0,0,1,0),
  Rating   = c(2.5,2,2,2.5,3,2,3,2.5,2,2.5,3,2.5)
)

res <- bm_gini_bm(df$Rating, df$Approval)
res
# $overall_gini ≈ 0.0833333
# $within_bm    ≈ 0.0282486
# $between_bm   ≈ 0.0550847
# $overlap      ≈ ~0  (groups do not overlap here)
# $check_sum    ≈ 0.0833333


df.res <- dineq::gini_decomp(x = df$Rating, z = df$Approval)

