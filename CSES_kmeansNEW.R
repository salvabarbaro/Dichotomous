kmax <- 4

fv.kmeans.per.id <- function(case_id, df, kmax.value = kmax) {
  
  df_case <- df %>%
    filter(case_ID == case_id)
  
  ids <- unique(df_case$id)
  
  optk_per_id <- sapply(ids, function(pid) {
    
    ratings <- df_case %>%
      filter(
        id == pid,
        !is.na(Rating),
        is.finite(Rating)
      ) %>%
      pull(Rating)
    
    if (
      length(ratings) < 6 ||
      length(unique(ratings)) < 2 ||
      !is.finite(sd(ratings)) ||
      sd(ratings) == 0
    ) {
      return(NA_integer_)
    }
    
    sc.rat <- as.matrix(scale(ratings))
    
    silh.data <- tryCatch(
      {
        factoextra::fviz_nbclust(
          x = sc.rat,
          FUNcluster = stats::kmeans,
          method = "silhouette",
          k.max = kmax.value,
          nstart = 25
        )[["data"]]
      },
      error = function(e) {
        return(NULL)
      }
    )
    
    if (
      is.null(silh.data) ||
      nrow(silh.data) == 0 ||
      all(is.na(silh.data$y))
    ) {
      return(NA_integer_)
    }
    
    # clusters is returned as a factor by factoextra
    silh.data$clusters <- as.integer(
      as.character(silh.data$clusters)
    )
    
    best_row <- which.max(silh.data$y)
    
    return(silh.data$clusters[best_row])
  })
  
  optk_clean <- optk_per_id[!is.na(optk_per_id)]
  
  table_k <- table(
    factor(
      optk_clean,
      levels = 2:kmax.value
    )
  )
  
  df_k <- as.data.frame(table_k)
  names(df_k) <- c("k", "count")
  
  df_k$k <- as.integer(as.character(df_k$k))
  df_k$case_ID <- case_id
  
  return(df_k)
}

case_ids <- unique(cs_p.df$case_ID)

# Use lapply to get list of data frames
res.list <- mclapply(
  case_ids,
  function(case_id) {
    fv.kmeans.per.id(
      case_id,
      cs_p.df,
      kmax.value = kmax
    )
  },
  mc.cores = 12
)

# Bind into one data frame
#new_kmeans <- dplyr::bind_rows(res.list)

new_kmeans <- bind_rows(res.list) %>%
  filter(k %in% c(2, 3, 4)) %>%
  pivot_wider(
    names_from = k,
    names_prefix = "k_",
    values_from = count,
    values_fill = 0
  ) %>%
  relocate(case_ID, .before = everything()) %>%
  mutate(
    total_k234 = k_2 + k_3 + k_4,
    k_2_pct = round(100 * k_2 / total_k234, 1),
    k_3_pct = round(100 * k_3 / total_k234, 1),
    k_4_pct = round(100 * k_4 / total_k234, 1)
  )


saveRDS(new_kmeans, "DATA/kmeansNEW.RDS")

## compare k-means(old) with k-means (new)
kmeans.res.old <- read.csv("DATA/csesSummary.csv", header = T)
kmeans.res.new <- readRDS("DATA/kmeansNEW.RDS")
k_compare.df <- kmeans.res.old %>%
  left_join(x = ., y = kmeans.res.new, by = "case_ID") %>%
  setNames(c("case_ID", 
"k.k_2", "k.k_3", "k.k_4", "k.total", "k.2pct", "k.3pct", "k.4pct",
"n.k_2", "n.k_3", "n.k_4", "n.total", "n.2pct", "n.3pct", "n.4pct"))
cor.test(k_compare.df$k.2pct, k_compare.df$n.2pct)

new_pam_compare.df <- kmeans.res.new %>%
    left_join(x = ., y = pam.res, by = "case_ID") %>%
  setNames(c("case_ID", 
"k.k_2", "k.k_3", "k.k_4", "k.total", "k.2pct", "k.3pct", "k.4pct",
"p.k_2", "p.k_3", "p.k_4", "p.total", "p.2pct", "p.3pct", "p.4pct"))
cor.test(new_pam_compare.df$k.2pct, new_pam_compare.df$p.2pct)



### Examples for the Appendix
cou.selection <- c("France_2012", "Germany_2021", "Great Britain_2019", "Israel_2020", "Japan_2017", "Tunisia_2019")

cou.table <- kmeans.res.new %>% filter(., case_ID %in% cou.selection) %>%
  dplyr::select(., c("case_ID", "k_2_pct", "k_3_pct", "k_4_pct")) 
stargazer::stargazer(cou.table, summary = F)
