#
library(ggplot2)
library(parallel)
library(scales)
library(latex2exp)
library(haven)
library(tidyverse)
library(modelsummary)
library(fixest)

cses <- read_dta("DATA/cses_imd.dta")
##
cses.df <- cses %>%
  mutate(
    C.Age = IMD2001_1,
    C.Gender = IMD2002,
    C.Education = IMD2003,
    C.Income = IMD2006,
    C.Ideology = IMD3006,
    C.SatDem = IMD3010,
    ID = IMD1005,
    Age = ifelse(C.Age > 99, NA, as.numeric(C.Age)),
    Gender = factor(case_when(
      C.Gender == 1 ~ "M",
      C.Gender == 2 ~ "F",
      TRUE ~ NA_character_
    )),
    Education = ifelse(C.Education > 4, NA, as.factor(C.Education)),
    IncomeQ = ifelse(C.Income > 5, NA, as.factor(C.Income)),
    Ideology = ifelse(C.Ideology > 10, NA, as.numeric(C.Ideology)), 
    SatisfactionDem = case_when(
      C.SatDem == 5 ~ 1,   # not at all satisfied
      C.SatDem == 4 ~ 2,   # not very satisfied
      C.SatDem == 6 ~ 3,   # neither nor
      C.SatDem == 2 ~ 4,   # fairly satisfied
      C.SatDem == 1 ~ 5,   # very satisfied
      TRUE ~ NA_real_
    ),   
    Country = IMD1006_NAM,
    Year = as.character(IMD1008_YEAR),
    case_ID = paste(Country, Year, sep = "_")
  ) %>%
  mutate(
    across(
      starts_with("IMD3008_"),
      ~ ifelse(.x < 11, .x, NA),
      .names = "party_rating_{tolower(sub('IMD3008_', '', .col))}"
    )) %>%
  mutate(Dist0 = abs(Ideology - 5)) %>%
  mutate(DistSq = Dist0^2) 
###########
rm(cses)  #
gc()      #
###########

optk_df <- readRDS("DATA/optkdf.RDS")   # generated through CSES_Cluster.R

# join back to original data
cses.all <- cses.df %>%
  left_join(optk_df, by = c("ID", "case_ID")) %>% 
  mutate(bin.k2 = ifelse(opt_k == 2, 1, 0))   # bin.k2 : LHS or the regressions

#### cses.all: the dataset for the regressions
###########################################################################################
### models with Dist0 as main variable, 
mod01 <- "bin.k2 ~ Dist0 + Age + Gender + Education + IncomeQ | case_ID"
mod02 <- "bin.k2 ~ Dist0 + Age + Gender | case_ID"
mod03 <- "bin.k2 ~ Dist0 + Age + IncomeQ | case_ID"
mod04 <- "bin.k2 ~ Dist0 + Education + IncomeQ   | case_ID"
mod05 <- "bin.k2 ~ Dist0  | case_ID"
#############################################
## Robustness Check: We replace Dist0 with SatisfactionDem in mod01
mod01.sat <- "bin.k2 ~ Satisfaction.Dem + Age + Gender + Education + IncomeQ | case_ID"
mod01.both <- "bin.k2 ~ Dist0 + Satisfaction.Dem + Age + Gender + Education + IncomeQ | case_ID"


feDist.fun <- function(m){
  feglm(
 fml = as.formula(m),
  family = binomial(link = "logit"),
  data   = cses.all
)
}

#distmod <- list(modD6, modD3, modD5, modD4, modD2, modD1)
distmod2 <- list(
  "Main" = mod01, 
  "Alt.01" = mod02, 
  "Alt.02" = mod03, 
  "Alt.03" = mod04, 
  "Alt.04" = mod05, 
  "Satisf." = mod01.sat, 
  "Both" = mod01.both)

DistReg <- lapply(distmod2, feDist.fun)

# robustness check

modelsummary(DistReg, 
  exponentiate = T, 
  stars = T, 
  statistic = "[{conf.low}, {conf.high}]",
  gof_map   = c("nobs", "aic", "bic"),
  vcov = "HC1",
  conf_level = 0.995,
    coef_rename = c(
    "Dist0" = "Ideol. Distance", 
    "IncomeQ" = "Income",
    "GenderM" = "Male",
    "Satisfaction.Dem" = "Satisf. Democ.")
#  output = "latex",
#  booktabs = TRUE,
#  file = "~/Documents/Research/Dichotomous/git/67b5f34c104b85acf4a11317/csesDistReg.tex"
)

modelplot(
  DistReg, 
  exponentiate = T, 
  vcov = "HC1", 
  conf_level = 0.995,
  coef_rename = c(
    "Dist0" = "Ideol. Distance", 
    "IncomeQ" = "Income",
    "GenderM" = "Male",
    "Satisfaction.Dem" = "Satisf. Democ.")
) + theme_bw(base_size = 24) + geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_viridis_d()
ggsave("RegressionPlots.pdf", width = 16, height = 8)


## relationship between Satisfaction.Dem and Dist0
ggplot(data = cses.all, aes(x = Satisfaction.Dem, y = Dist0))  +
  geom_hex() +
  theme_bw(base_size = 22)

car::vif(lm(bin.k2 ~ Satisfaction.Dem + Dist0, data = cses.all))
## close to 1: no multicollinearity problem

prop.table(table(cses.all$Satisfaction.Dem, cses.all$Dist0))








### Regression analysis
mod01 <- "bin.k2 ~ Age + Gender + Education + IncomeQ + Ideology + Satisfaction.Dem | case_ID"
mod02 <- "bin.k2 ~ Age + Gender + Education  | case_ID"
mod03 <- "bin.k2 ~ Education + Ideology  | case_ID"
mod04 <- "bin.k2 ~ Education + Ideology + Satisfaction.Dem  | case_ID"
mod05 <- "bin.k2 ~ Education + Satisfaction.Dem  | case_ID"
mod06 <- "bin.k2 ~ Education  | case_ID"

### models with Democracy Satisfaction as main variable
mod01 <- "bin.k2 ~ Age + Gender + Education + IncomeQ + Ideology + Satisfaction.Dem | case_ID"
mod02 <- "bin.k2 ~ Age + Gender + Satisfaction.Dem  | case_ID"
mod03 <- "bin.k2 ~ Satisfaction.Dem + Ideology  | case_ID"
mod04 <- "bin.k2 ~ Education + Ideology + Satisfaction.Dem  | case_ID"
mod05 <- "bin.k2 ~ Education + Satisfaction.Dem  | case_ID"
mod06 <- "bin.k2 ~ Satisfaction.Dem  | case_ID"



felogreg.fun <- function(m){
  feglm(
 fml = as.formula(m),
  family = binomial(link = "logit"),
  data   = cses.df
)
}

mods <- list(mod01, mod02, mod03, mod04, mod05, mod06)
mods <- list(mod06, mod05, mod03, mod02, mod04, mod01)
logregs <- lapply(mods, felogreg.fun)

modelsummary(logregs, 
  exponentiate = T, 
  stars = T, 
  statistic = '[{conf.low}, {conf.high}]',
  gof_omit = "R2|RMSE",
#  gof_map   = c("nobs", "aic", "bic"),
  vcov = "HC2",
  conf_level = 0.95)

options("modelsummary_format_numeric_latex" = "plain")
modelsummary(
  logregs, 
  exponentiate = TRUE, 
  stars = TRUE, 
  statistic = "[{conf.low}, {conf.high}]",
  gof_map   = c("nobs", "aic", "bic"),
  vcov = "HC1",
  conf_level = 0.95,
  output = "latex",
  booktabs = TRUE,
  file = "~/Documents/Research/Dichotomous/git/67b5f34c104b85acf4a11317/csesLogReg.tex"
)

#### Dichotomous preferences and Ideology.
cses.sat <- cses.df %>% dplyr::select(., c(bin.k2, Ideology, Satisfaction.Dem))
ideo.df <- cses.sat %>% group_by(Ideology) %>% reframe(DP =mean(bin.k2, na.rm = T))
sati.df <- cses.sat %>% group_by(Satisfaction.Dem) %>% reframe(DP =mean(bin.k2, na.rm = T))
cses.sat <- cses.df %>%
  rowwise() %>%
  mutate(
    count_0_10 = sum(c_across(starts_with("party_rating")) %in% c(0, 10), na.rm = TRUE),
    count_NA   = sum(is.na(c_across(starts_with("party_rating"))))
  ) %>%
  ungroup()
cses.sat %>% group_by(Satisfaction.Dem) %>%
  reframe(ext = mean(count_0_10, na.rm = T),
          na  = mean(count_NA,   na.rm = T)
          )

write_rds(cses.sat, file = "DATA/csesDemSat.RDS")
cses.sat <- readRDS("DATA/csesDemSat.RDS")

glm.ideology <- feglm(fml = bin.k2 ~ Ideology, family = binomial(link = "logit"),
  data   = cses.df)
modelsummary(glm.ideology, 
  exponentiate = T, 
  stars = T, 
  statistic = "[{conf.low}, {conf.high}]",
  gof_map   = c("nobs", "aic", "bic"),
  vcov = "HC1",
  conf_level = 0.95)

cses.dist <- cses.df %>%
  mutate(Dist0 = abs(Ideology - 5)) %>%
  mutate(DistSq = Dist0^2) %>%
  mutate(bin.k2 = ifelse(opt_k == 2, 1, 0))

modD1 <- "bin.k2 ~ Age + Gender + Education + IncomeQ + Dist0 + Satisfaction.Dem | case_ID"
modD2 <- "bin.k2 ~ Age + Gender + Education + DistSq  | case_ID"
modD3 <- "bin.k2 ~ Education + Dist0  | case_ID"
modD4 <- "bin.k2 ~ Education + Dist0 + Satisfaction.Dem  | case_ID"
modD5 <- "bin.k2 ~ Education + Satisfaction.Dem + DistSq  | case_ID"
modD6 <- "bin.k2 ~ Education + Dist0  | case_ID"





cses.dist %>% group_by(Satisfaction.Dem) %>%
  reframe(mnDist = mean(Dist0, na.rm = T),
          mn.DSq = mean(DistSq, na.rm = T))



## Dissatisfaction and Idelo Distance may be too correlated to each other
ggplot(
  cses.dist,
  aes(
    x = Ideology,
    y = Satisfaction.Dem
  )
) +
  geom_jitter(
    alpha = 0.2,
    width = 0.1,
    height = 0.1
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "red"
  ) +
  theme_bw(base_size = 18)

#cs_p_wide <- cs_p.df %>%
#  pivot_wider(
#    id_cols = c(case_ID, ID),      # keep respondent + election IDs
#    names_from = party,            # variable that held party names
#    values_from = Rating,          # the ratings
#    names_prefix = "party_rating_" # optional, keep consistent naming
#  )
#colnames(cs_p_wide)









## Cluster analyses on CSES data
# 1. Load data
#load("~/Documents/Research/Elections/AnnaProjects/CondorcetParadox/Data/cses_imd.rdata")
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
## Regressions on variables
## IMD5013       >>> ELECTORAL FORMULA IN ALL SEGMENTS: LOWER HOUSE
##  1. MAJORITARIAN  2. PROPORTIONAL  3. MIXED  9. MISSING
cses.all <- cses.all %>% 
  mutate(
    C.ElectSystem = ifelse(IMD5013 == 9, NA, as.factor(IMD5013)),
    NbEffParties  = ifelse(IMD5058_1 > 100, NA, as.numeric(IMD5058_1)) ) 

cses.short <- cses.all %>% 
  dplyr::select(., c("case_ID", "C.ElectSystem", "NbEffParties", "IMD1006_UNALPHA3")) %>% unique(.)

cses.elecsystem <- summary_tables %>%
  left_join(
    x = ., 
    y = cses.short,
    by = "case_ID"
  )

reg.es <- lm(
  formula = k_2_pct ~ as.factor(C.ElectSystem) + NbEffParties, 
  data = cses.elecsystem
)
summary(reg.es)

library(estimatr)
reg.esClst <- lm_robust(
  formula = k_2_pct ~ as.factor(C.ElectSystem) + NbEffParties, 
  data = cses.elecsystem, 
  clusters = IMD1006_UNALPHA3
)
summary(reg.esClst)





























