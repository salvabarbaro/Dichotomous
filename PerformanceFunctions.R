#### Replication file to
#### Dichotomous Preferences: Concepts, Measurement, and Evidence
#### by Salvatore Barbaro and Anna-Sophie Kurella
#### This replication file covers all analyses with the CSES dataset
#### Section 5
######################################################################################
library(tidyverse)
library(cluster)
library(factoextra)
library(tidyr)
#library(parallel)
library(future)
library(future.apply)
library(scales)
library(latex2exp)
library(fixest)
library(modelsummary)
#######################################################################################
## 1.  Read and prepare data
#######################################################################################
# Get the original CSES-data and recode/rename variables of interest
cses_imd <- haven::read_dta("DATA/cses_imd.dta")
cses.df <- cses_imd %>%
  mutate(C.ID = IMD1005,
         C.Age = IMD2001_1,
         C.Gender = IMD2002,
         C.Education = IMD2003,
         C.Income = IMD2006,
         C.SocEconStatus = IMD2016,
         C.Ideology = IMD3006,
         C.SatisfDem = IMD3010,
         C.SatDem = IMD3010,    # the same as C.SatisfDem
         C.ElectSystem = IMD5013,
         C.NbEffParties = IMD5058_1,
         C.case_ID = paste(IMD1006_NAM, IMD1008_YEAR, sep = "_")
        ) %>% 
  dplyr::select(starts_with("C.")) %>%
  mutate(
    ID = C.ID,
    Age = ifelse(C.Age > 99, NA, as.numeric(C.Age)),
    Gender = factor(case_when(
      C.Gender == 1 ~ "M",
      C.Gender == 2 ~ "F",
      TRUE ~ NA_character_  )),
    Education = ifelse(C.Education > 4, NA, as.factor(C.Education) ),
    Income = ifelse(C.Income > 5, NA, as.factor(C.Income)),
    IncomeQ = ifelse(C.Income > 5, NA, as.factor(C.Income)),
    Ideology =  ifelse(C.Ideology > 10, NA, C.Ideology),
    Dissatisfaction = ifelse(C.SatisfDem >5, NA, as.factor(C.SatisfDem)),
    SatisfactionDem = case_when(
      C.SatDem == 5 ~ 1,   # not at all satisfied
      C.SatDem == 4 ~ 2,   # not very satisfied
      C.SatDem == 6 ~ 3,   # neither nor
      C.SatDem == 2 ~ 4,   # fairly satisfied
      C.SatDem == 1 ~ 5,   # very satisfied
      TRUE ~ NA_real_
    ),
    ElectSystem = ifelse(C.ElectSystem == 9, NA, as.factor(C.ElectSystem)),
    NbEffParties  = ifelse(C.NbEffParties > 100, NA, as.numeric(C.NbEffParties)),
    case_ID = C.case_ID,
    Country = cses_imd$IMD1006_UNALPHA3,
    Year = as.character(cses_imd$IMD1008_YEAR),
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

## clean-up
rm(cses_imd)
gc()
#
# The second cses-data is taken from Barbaro/Kurella: Condorcet Paradox (Public Choice)
## It effectively uses the same cses_imd as above and renames the variables regarding
## party ratings. Please find cses at the github-repo 
## https://github.com/salvabarbaro/CondorcetCycle/tree/main/Data
#https://cses.org/wp-content/uploads/2024/02/cses_imd_codebook_part2_variables.txt
cses <- readRDS("DATA/cses.RDS")
## Approach: we consider only individuals who rated at least six parties.
## Given this restriction, we consider only elections with at least 100 individuals
## This effectively removes case_IDs with only few parties considered. 
## Effects: 212.729 participants, 172 elections. 
cs_p.df <- cses %>% 
  pivot_longer(cols = starts_with("party_rating")) %>% 
  dplyr::rename(Rating = value) %>%
  group_by(case_ID, id) %>%
  filter(sum(!is.na(Rating)) >= 6) %>%
  ungroup() %>%
  group_by(case_ID) %>%
  filter(n_distinct(id) >= 100) %>%
  ungroup()

# clean-up
rm(cses)
gc()

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
################################################################################
## 2. Cluster Analyses 
################################################################################
## we will use mclapply. This function works (exclusively) on Linux. For replication
## with other OS (like win, iOS) you need to adjust the parallelisation part or 
## remove it. Serious damage can happen by running mclapply() on non-linux OS!
## To avoid parallelisation, replace mclapply() by lapply() (Takes about 8 hours to run)
#
## Define the number of cores you want to consider for parallelisation
nb.cores <- 18
## Define the max number of clusters 
k_max <- 4
## Define the number of iterations (starting points).
nb.iter <- 25  # 25 ist a common number in research
##
set.seed(55234)
##

kmeans.id <- function(df_id, kmax = k_max, nstart = nb_iter) {

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

    km <- stats::kmeans(
      x,
      centers = k,
      nstart = nstart
    )

    sil <- mean(
      cluster::silhouette(
        km$cluster,
        d
      )[, "sil_width"]
    )

    list(
      k = k,
      fit = km,
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

  df_id %>%
    mutate(
      opt_k = opt_k,
      cluster = final_fit$cluster,
      cluster_center =
        final_fit$centers[final_fit$cluster, 1]
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

future::plan(
  future::multisession,
  workers = nb.cores
)

res <- future.apply::future_lapply(
  split.df,
  kmeans.id,
  kmax = 4,
  nstart = 25,
  future.seed = TRUE
)

#res <- parallel::mclapply(
#  split.df,
#  kmeans.id,
#  kmax = 4,
#  nstart = 25,
#  mc.cores = nb.cores
#)

cluster.df <- dplyr::bind_rows(res)
cluster.short <- cluster.df %>% dplyr::select(., c("id", "opt_k", "case_ID")) %>% distinct()
# test of consistency:
cluster.short %>%
  count(id) %>%
  count(n)
# test of consistency - end 

cou.table <- cluster.short %>%
  mutate(
    opt_k_cat = ifelse(opt_k >= 4, 4, opt_k)
  ) %>%
  count(case_ID, opt_k_cat) %>%
  tidyr::pivot_wider(
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

cou.selection <- c("France_2012", "Germany_2021", "Great Britain_2019", "Israel_2020", "Japan_2017", "Tunisia_2019")

## some examples
cou.appendixtable <- cou.table %>% dplyr::filter(., case_ID %in% cou.selection) %>%
  dplyr::select(., c("case_ID", "k2_pct", "k3_pct", "k4_pct")) %>%
  tidyr::extract(
    case_ID,
    into = c("Country", "Year"),
    regex = "^(.*)_([0-9]{4})$",
    remove = FALSE
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  dplyr::select(., -c("case_ID")) %>%
  setNames(c("Country", "Year", "k=2", "k=3", "k=4"))
#
# all countries
cou.fulltable <- cou.table %>% 
  dplyr::select(., c("case_ID", "k2_pct", "k3_pct", "k4_pct")) %>%
  tidyr::extract(
    case_ID,
    into = c("Country", "Year"),
    regex = "^(.*)_([0-9]{4})$",
    remove = FALSE
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  dplyr::select(., -c("case_ID")) %>%
  setNames(c("Country", "Year", "k=2", "k=3", "k=4"))


kableExtra::kbl(
  cou.appendixtable,
  format = "latex",
  booktabs = TRUE,
  digits = 2,
  label = "tb.csesexamples",
  caption = "Results from some selected countries/elections", 
  linesep = NULL
)

# for the appendix: long table
kableExtra::kbl(
  cou.fulltable,
  format = "latex",
  booktabs = TRUE,
  digits = 2,
  label = "tb.csesexamples",
  caption = "Results from some selected countries/elections", 
  linesep = NULL
)

res_long <- cou.table %>%
  pivot_longer(cols = ends_with("pct"),
               names_to = "k_type",
               values_to = "percent") %>%
  mutate(k_type = recode(k_type,
                         "k2_pct" = "k2",
                         "k3_pct" = "k3",
                         "k4_pct" = "k4"))


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
  theme_minimal(base_size = 24) +
  theme(legend.position = "none")
ggsave("cses169.pdf", plot = p1, width = 16, height = 8)
rm(p1, res_long)
gc()
#################################################################
## 3. Regression
##################################################################
### Regression: Can we explain the share of k2?
## Polarization data
polariz.df <- readRDS(file = "polarization.RDS") 
# adding polarization data to cses.df
cses.df <- cses.df %>% 
  rename("caseIDiso3c" = "case_ID") %>%
  rename("case_ID" = "C.case_ID")   %>%
  left_join(x = ., y = polariz.df, by = "case_ID") %>%   
  left_join(x = ., y = cou.table, by = "case_ID")

cses.ols <- cses.df %>% 
  dplyr::select(., c("k2_pct", "ElectSystem", "NbEffParties", 
  "case_ID", "Country",
  "polarization_parties", "polarization_voter")) %>%
  unique() %>%
  mutate(year = as.integer(stringr::str_extract(case_ID, "\\d{4}$")))

regmods.k2 <- list(
  "Main" = "k2_pct ~ polarization_parties + factor(ElectSystem) + NbEffParties + year",
  "RC1"   = "k2_pct ~  factor(ElectSystem) + NbEffParties + year",
  "RC2" = "k2_pct ~ polarization_parties + factor(ElectSystem) + NbEffParties"
)

olsfun <- function(m){
  reg = lm(formula = m, data = cses.ols)
}

olsregs <- lapply(regmods.k2, olsfun)

options("modelsummary_format_numeric_latex" = "plain")
modelsummary::modelsummary(olsregs, 
  stars = T, 
#  gof_omit = "AIC|BIC|Log.|RMSE",
  coef_map = c(
    "polarization_parties" = "Party Polarization",
    "polarization_voter" = "Voter Polarization",
    "NbEffParties" = "Nb. Eff. Parties",
    "factor(ElectSystem)2" = "ElecSys-Proport.",
    "factor(ElectSystem)3" = "ElecSys-Mixed",
    "year" = "Time",
    "Country" = "Country" ),
  vcov = list(
    ~ Country,
    ~ Country,
    ~ Country + year),
  gof_omit = "AIC|BIC|Log.|RMSE"#, 
#  output = "k2ols.tex"
)

cses.cre <- cses.ols %>%
  group_by(Country) %>%
  mutate(
    polarization_between = mean(polarization_parties, na.rm = TRUE),
    polarization_within =
      polarization_parties - polarization_between,

    parties_between = mean(NbEffParties, na.rm = TRUE),
    parties_within =
      NbEffParties - parties_between
  ) %>%
  ungroup()
#
m.cre <- lme4::lmer(
  k2_pct ~ polarization_within +
    polarization_between +
    factor(ElectSystem) +
    parties_within +
    parties_between +
    year +
    (1 | Country),
  data = cses.cre
)

gof_map <- tribble(
  ~raw,          ~clean,                 ~fmt,
  "nobs",        "Observations",         0,
  "ngrps",       "Countries",            0,
  "var__Country","Var(Random intercept)",2,
  "var__Residual","Var(Residual)",       2
)

coef_map = c(
    "polarization_within" = "Polarization (Within)",
    "polarization_between" = "Polarization (between)",
    "parties_within" = "Nb. Eff. Parties (within)",
    "parties_between" = "Nb. Eff. Parties (between)",
    "factor(ElectSystem)2" = "ElecSys-Proport.",
    "factor(ElectSystem)3" = "ElecSys-Mixed",
    "year" = "Time",
    "Country" = "Country" )

modelsummary(
  m.cre, stars = T,
  gof_map = gof_map,
  coef_map = coef_map#,
#  output = "Mundlak.tex"
)
################################################
# Micro-Level
optk_df <- cluster.short %>% rename(ID = id)
cses.all <- cses.df  %>%
  left_join(x = ., y = optk_df, by = c("case_ID", "ID")) %>% 
  mutate(bin.k2 = ifelse(opt_k == 2, 1, 0))   # bin.k2 : LHS or the regressions

#
### models with Dist0 as main variable, 
mod01 <- "bin.k2 ~ Dist0 + Age + Gender + Education + IncomeQ | case_ID"
mod02 <- "bin.k2 ~ Dist0 + Age + Gender | case_ID"
mod03 <- "bin.k2 ~ Dist0 + Age + IncomeQ | case_ID"
mod04 <- "bin.k2 ~ Dist0 + Education + IncomeQ   | case_ID"
mod05 <- "bin.k2 ~ Dist0  | case_ID"
#############################################
## Robustness Check: We replace Dist0 with SatisfactionDem in mod01
mod01.sat <- "bin.k2 ~ SatisfactionDem + Age + Gender + Education + IncomeQ | case_ID"
mod01.both <- "bin.k2 ~ Dist0 + SatisfactionDem + Age + Gender + Education + IncomeQ | case_ID"


feDist.fun <- function(m){
  fixest::feglm(
 fml = as.formula(m),
  family = binomial(link = "logit"),
  data   = cses.all  )}

micromodels <- list(
  "Main" = mod01, 
  "Alt.01" = mod02, 
  "Alt.02" = mod03, 
  "Alt.03" = mod04, 
  "Alt.04" = mod05, 
  "Satisf." = mod01.sat, 
  "Both" = mod01.both)

DistReg <- lapply(micromodels, feDist.fun)

coef_map = c(
    "SatisfactionDem" = "Satisf. Democracy",
    "IncomeQ" = "Income",
    "Education" = "Education",
    "GenderM" = "Gender (Male)",
    "Age" = "Age",
    "Dist0" = "Ideology")


modelsummary(
  DistReg, 
  exponentiate = T,
  stars = T,
  conf_level = 0.9,
  statistic = '[{conf.low}, {conf.high}]',
  gof_omit = "R2|RMSE",
  vcov = ~case_ID, 
  coef_map = coef_map)

modelplot(
  DistReg, 
  exponentiate = T,
  vcov = ~case_ID,
  coef_map = coef_map,
  conf_level = 0.99
) + theme_bw(base_size = 22) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_colour_viridis_d()
ggsave("microregressions.pdf", width = 16, height = 9)
###################################################################################################################################