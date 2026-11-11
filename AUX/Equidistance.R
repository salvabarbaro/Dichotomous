## Checking for equidistance
## requires: cs_p.df from RF_CSES.R
library(tidyverse)

equidistance.df <- cs_p.df %>%
  group_by(case_ID, id) %>%
  summarise(
    n_candidates = sum(!is.na(Rating)),
    ratings = list(sort(Rating[!is.na(Rating)])),
    .groups = "drop"
  ) %>%
  mutate(
    exact_equidistant = sapply(ratings, function(x) {
      if (length(x) < 3) return(FALSE)

      d <- diff(x)

      # All distances must be identical and strictly positive
      all(d == d[1]) && d[1] > 0
    })
  )

# Overall frequency
equidistance.df %>%
  summarise(
    N = n(),
    n_equidistant = sum(exact_equidistant),
    pct_equidistant = 100 * mean(exact_equidistant)
  )

# Frequency by election
equidistance.df %>%
  group_by(case_ID) %>%
  summarise(
    N = n(),
    n_equidistant = sum(exact_equidistant),
    pct_equidistant = 100 * mean(exact_equidistant),
    .groups = "drop"
  )
