## CSES-data with PAM
library(dplyr)
library(factoextra)
library(cluster)
library(parallel)

kmax <- 4

fv.pam.per.id <- function(case_id, df, kmax.value = kmax) {
  
  df_case <- df %>%
    filter(case_ID == case_id)
  
  ids <- unique(df_case$id)
  
  optk_per_id <- sapply(ids, function(pid) {
    
    ratings <- df_case %>%
      filter(id == pid, !is.na(Rating)) %>%
      pull(Rating)
    
    # Mindestzahl verfügbarer Ratings
    if (length(ratings) < 6) {
      return(NA_integer_)
    }
    
    # Standardisierung
    sc.rat <- scale(ratings)
    sc.rat <- as.matrix(sc.rat)
    
    silh.data <- tryCatch(
      {
        factoextra::fviz_nbclust(
          x = sc.rat,
          FUNcluster = cluster::pam,
          method = "silhouette",
          k.max = kmax.value
        )[["data"]]
      },
      error = function(e) {
        return(NULL)
      }
    )
    
    if (is.null(silh.data) || all(is.na(silh.data$y))) {
      return(NA_integer_)
    }
    
    # k mit der höchsten durchschnittlichen Silhouettenbreite
    best_row <- which.max(silh.data$y)
    
    return(as.integer(silh.data$clusters[best_row]))
  })
  
  optk_clean <- na.omit(optk_per_id)
  
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
res.list <- mclapply(case_ids, function(case_id) {
  fv.pam.per.id(case_id, cs_p.df, kmax.value = kmax)
}, mc.cores = 18)

# Bind into one data frame
summary_tables <- dplyr::bind_rows(res.list)

summary_tables <- bind_rows(res.list) %>%
  filter(k %in% c("2", "3", "4")) %>%      # only keep k = 2, 3, 4
  pivot_wider(
    names_from = k,
    names_prefix = "k_",
    values_from = count,
    values_fill = 0                      # fill missing with 0
  ) %>%
  relocate(case_ID, .before = everything())  %>% # case_ID as first column
  mutate(
    total_k234 = k_2 + k_3 + k_4,
    k_2_pct = round(100 * k_2 / total_k234, 1),
    k_3_pct = round(100 * k_3 / total_k234, 1),
    k_4_pct = round(100 * k_4 / total_k234, 1)
  )

head(summary_tables)
saveRDS(summary_tables, "DATA/PAMSummary.rds")
pam.res <- readRDS("DATA/PAMSummary.rds")


## compare with k-means clustering
kmeans.res <- read.csv("DATA/csesSummary.csv", header = T)

common.df <- kmeans.res %>%
  left_join(x = ., y = pam.res, by = "case_ID") %>%
  setNames(c("case_ID", 
"k.k_2", "k.k_3", "k.k_4", "k.total", "k.2pct", "k.3pct", "k.4pct",
"p.k_2", "p.k_3", "p.k_4", "p.total", "p.2pct", "p.3pct", "p.4pct"))

colnames(common.df)
cor(common.df$k.2pct, common.df$p.2pct)
ggplot(data = common.df, aes(x = k.2pct, y = p.2pct)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_gray(base_size = 22) +
  labs(x = "k-means", y = "PAM")

cor.test(common.df$k.2pct, common.df$p.2pct)
