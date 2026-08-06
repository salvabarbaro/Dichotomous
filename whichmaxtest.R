library(dplyr)
library(factoextra)

kmax <- 4
case_id <- "Germany_2021"

eligible_ids <- cs_p.df %>%
  filter(
    case_ID == case_id,
    !is.na(Rating),
    is.finite(Rating)
  ) %>%
  group_by(id) %>%
  summarise(
    n_ratings = n(),
    n_distinct_ratings = n_distinct(Rating),
    rating_sd = sd(Rating),
    .groups = "drop"
  ) %>%
  filter(
    n_ratings >= 6,
    n_distinct_ratings >= kmax,
    is.finite(rating_sd),
    rating_sd > 0
  )

eligible_ids

pid <- eligible_ids$id[1]

ratings <- cs_p.df %>%
  filter(
    case_ID == case_id,
    id == pid,
    !is.na(Rating),
    is.finite(Rating)
  ) %>%
  pull(Rating)

ratings

length(ratings)
length(unique(ratings))
sd(ratings)
anyNA(ratings)

sc.rat <- scale(ratings)
sc.rat <- as.matrix(sc.rat)

anyNA(sc.rat)
any(!is.finite(sc.rat))

set.seed(55234)

tmp <- factoextra::fviz_nbclust(
  sc.rat,
  kmeans,
  method = "silhouette",
  k.max = kmax
)

tmp$data

silh.values <- tmp$data$y

old_result <- which.max(silh.values)

old_result


###########################################################
check_k_values <- function(case_id, df, kmax.value = 4) {

  df_case <- df %>%
    filter(case_ID == case_id)

  ids <- unique(df_case$id)

  results <- lapply(ids, function(pid) {

    ratings <- df_case %>%
      filter(
        id == pid,
        !is.na(Rating),
        is.finite(Rating)
      ) %>%
      pull(Rating)

    if (length(ratings) < 6 || sd(ratings) == 0) {
      return(NULL)
    }

    sc.rat <- as.matrix(scale(ratings))

    silh.data <- tryCatch(
      factoextra::fviz_nbclust(
        x = sc.rat,
        FUNcluster = stats::kmeans,
        method = "silhouette",
        k.max = kmax.value,
        nstart = 25
      )$data,
      error = function(e) NULL
    )

    if (is.null(silh.data)) {
      return(
        data.frame(
          id = pid,
          error = TRUE
        )
      )
    }

    data.frame(
      id = pid,
      k = silh.data$clusters,
      silhouette = silh.data$y,
      error = FALSE
    )
  })

  bind_rows(results)
}

check.deu <- check_k_values(
  "Germany_2021",
  cs_p.df,
  kmax.value = 4
)
table(check.deu$k, useNA = "ifany")

optimal.deu <- check.deu %>%
  filter(
    !is.na(k),
    !is.na(silhouette),
    is.finite(silhouette)
  ) %>%
  group_by(id) %>%
  slice_max(
    order_by = silhouette,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

table(optimal.deu$k)
prop.table(table(optimal.deu$k))

fv.kmeans.per.id("Germany_2021", cs_p.df)

class(tmp$data$clusters)
levels(tmp$data$clusters)

best_row <- which.max(tmp$data$y)

tmp$data$clusters[best_row]
as.integer(tmp$data$clusters[best_row])
as.integer(as.character(tmp$data$clusters[best_row]))
