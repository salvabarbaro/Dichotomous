library(dplyr)
library(cluster)
library(factoextra)
library(ggplot2)
library(ineq)
library(rstatix)
library(tidyr)
library(IC2)
library(parallel)
library(scales)
library(latex2exp)
library(haven)

## Cluster analyses on CSES data
setwd("~/Documents/Research/Dichotomous/github/Dichotomous/DATA/")
# 1. Load data
cses_imd <- haven::read_dta("~/Documents/Research/Dichotomous/Data/cses_imd_stata/cses_imd.dta")
cses.df <- cses_imd %>%
  mutate(C.ID = IMD1005,
         C.Age = IMD2001_1,
         C.Gender = IMD2002,
         C.Education = IMD2003,
         C.Income = IMD2006,
         C.SocEconStatus = IMD2016,
         C.Ideology = IMD3006,
         C.SatisfDem = IMD3010) %>%
  dplyr::select(starts_with("C.")) %>%
  mutate(
    Age = ifelse(C.Age > 99, NA, as.numeric(C.Age)),
    Gender = factor(case_when(
      C.Gender == 1 ~ "M",
      C.Gender == 2 ~ "F",
      TRUE ~ NA_character_  )),
    Education = ifelse(C.Education > 4, NA, as.factor(C.Education) ),
    Income = ifelse(C.Income > 5, NA, as.factor(C.Income)),
    Ideology =  ifelse(C.Ideology > 10, NA, C.Ideology),
    Dissatisfaction = ifelse(C.SatisfDem >5, NA, as.factor(C.SatisfDem))
  )

#load("~/Documents/Research/Elections/AnnaProjects/CondorcetParadox/Data/cses_imd.rdata")
# Codebook:
#https://cses.org/wp-content/uploads/2024/02/cses_imd_codebook_part2_variables.txt
cses <- read.csv("cses.csv", header = T)
## Approach: we consider only individuals who rated at least six parties.
## Given this restriction, we consider only elections with at least 100 individuals
## This effectively removes case_IDs with only few parties considered. 
## Effects: 212.729 participants, 172 elections. 
cs_p.df <- cses %>% 
  dplyr::select(-starts_with("candidate_rating")) %>%
  pivot_longer(cols = starts_with("party_rating")) %>% 
  dplyr::rename(Rating = value) %>%
  group_by(case_ID, id) %>%
  filter(sum(!is.na(Rating)) >= 6) %>%
  ungroup() %>%
  group_by(case_ID) %>%
  filter(n_distinct(id) >= 100) %>%
  ungroup()
#
length(unique(cs_p.df$case_ID))  # number of elections
length(unique(cs_p.df$id))       # number of respondents
cs_p.df %>%
  mutate(Country = sub("_.*", "", case_ID)) %>%   # take everything before the underscore
  summarise(n_countries = n_distinct(Country))  # number of countries
respondents_summary <- cs_p.df %>%
  group_by(case_ID) %>%
  summarise(n_respondents = n_distinct(id), .groups = "drop")
respondents_summary %>%
  summarise(
    min     = min(n_respondents),
    mean    = mean(n_respondents),
    max     = max(n_respondents),
    n_cases = n()  # nb elections
  )
rm(respondents_summary)
#################################################################################################
## Cluster-analytical function
kmax <- 4
###############
fv.per.id <- function(case_id, df, kmax.value = kmax) {
  df_case <- df %>% filter(case_ID == case_id)
  ids <- unique(df_case$id)
  optk_per_id <- sapply(ids, function(pid) {
    ratings <- df_case %>% filter(id == pid, !is.na(Rating)) %>% pull(Rating)
    if(length(ratings) < 6) return(NA)
    eps <- sample(seq(-.1, .1, .001), 1)
    rat.temp <- ratings + eps
    sc.rat <- scale(rat.temp)
    sc.rat <- as.matrix(sc.rat)
    silh.values <- tryCatch({
      fviz_nbclust(sc.rat, kmeans, method = "silhouette", k.max = kmax)[["data"]][["y"]]
    }, error = function(e) {
      return(rep(NA, kmax.value))
    })

    if(all(is.na(silh.values))) return(NA)

    return(which.max(silh.values))
  })

  optk_clean <- na.omit(optk_per_id)

  # Tabulate only existing k values
  table_k <- table(factor(optk_clean, levels = 2:kmax.value))

  # Convert to data frame
  df_k <- as.data.frame(table_k)
  names(df_k) <- c("k", "count")
  df_k$case_ID <- case_id

  return(df_k)
}

####
case_ids <- unique(cs_p.df$case_ID)

# Use lapply to get list of data frames
res.list <- mclapply(case_ids, function(case_id) {
  fv.per.id(case_id, cs_p.df, kmax.value = kmax)
}, mc.cores = 8)

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

write.csv(summary_tables, "csesSummary.csv", row.names = F)

summary_tables <- read.csv("csesSummary.csv", header = T)
sapply(res[,6:8], mean)

res_long <- res %>%
  pivot_longer(cols = ends_with("pct"),
               names_to = "k_type",
               values_to = "percent") %>%
  mutate(k_type = recode(k_type,
                         "k_2_pct" = "k2",
                         "k_3_pct" = "k3",
                         "k_4_pct" = "k4"))


p1 <- ggplot(res_long, aes(x = k_type, y = percent)) +
  geom_boxplot(fill = "steelblue", colour = "darkblue", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, colour = "darkred") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(labels = c(
    "k2" = TeX("$\\tilde{k}=2$"),
    "k3" = TeX("$\\tilde{k}=3$"),
    "k4" = TeX("$\\tilde{k}=4$")
  )) +
  labs(x = NULL, y = "Percentage") +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none")
ggsave("~/Documents/Research/Dichotomous/git/67b5f34c104b85acf4a11317/cses.pdf", plot = p1, width = 16, height = 9)



###########################################################




























