#### Replication file to
#### Dichotomous Preferences...
#### by Barbaro and Kurella
#### This replication file covers all analyses with the CSES dataset
#### Section 5
######################################################################################
library(tidyverse)
library(cluster)
library(factoextra)
library(tidyr)
library(parallel)
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

#
# The second cses-data is taken from Barbaro/Kurella: Condorcet Paradox (Public Choice)
#https://cses.org/wp-content/uploads/2024/02/cses_imd_codebook_part2_variables.txt
cses <- read.csv("DATA/cses.csv", header = T)
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
################################################################################
## 2. Cluster Analyses 
################################################################################
## we will use mclapply. This function works (exclusively) on Linux. For replication
## with other OS (like win, iOS) you need to adjust the parallelisation part or 
## remove it. Serious damage can happen by running mclapply() on non-linux OS!
## To avoid parallelisation, replace mclapply() by lapply() (Takes about 8 hours to run)
#
## Define the number of cores you want to consider for parallelisation
nb.cores <- 12
## Define the max number of clusters 
kmax <- 4
## Define the number of iterations (starting points).
nb.iter <- 25  # 25 ist a common number in research
##

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
          nstart = nb.iter
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
  mc.cores = nb.cores
)

## the term new_kmeans refers to the fact that in the initial version we used a 
## smaller number of initial iterations (nstart) and adjusted the number in
## response to a reviewer comment

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
## for quick replications: read the kmeansNEW.RDS
new_kmeans <- readRDS("DATA/kmeansNEW.RDS")
#################################################################################
## 2.1 In the appendix, we displayed some examples to capture the distribution 
## of optimal cluster numbers in some elections
## table for case studies
cou.selection <- c("France_2012", "Germany_2021", "Great Britain_2019", "Israel_2020", "Japan_2017", "Tunisia_2019")

cou.table <- new_kmeans %>% filter(., case_ID %in% cou.selection) %>%
  dplyr::select(., c("case_ID", "k_2_pct", "k_3_pct", "k_4_pct")) %>%
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
  cou.table,
  format = "latex",
  booktabs = TRUE,
  digits = 2,
  label = "tb.csesexamples",
  caption = "Results from some selected countries/elections",
  linesep = NULL
)
rm(cou.selection, cou.table)
###########################################################################################
## Display results: Figure 2
res_long <- new_kmeans %>%
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
  left_join(x = ., y = new_kmeans, by = "case_ID")

cses.ols <- cses.df %>% 
  dplyr::select(., c("k_2_pct", "ElectSystem", "NbEffParties", 
  "case_ID", "Country",
  "polarization_parties", "polarization_voter")) %>%
  unique() %>%
  mutate(year = as.integer(stringr::str_extract(case_ID, "\\d{4}$")))

#saveRDS(cses.ols, file = "DATA/csesOLS.RDS")
#cses.ols <- readRDS("DATA/csesOLS.RDS") 
#### Regressions
#### Makro-Level
## models
regmods.k2 <- list(
  "Main" = "k_2_pct ~ polarization_parties + factor(ElectSystem) + NbEffParties + year",
  "RC1"   = "k_2_pct ~  factor(ElectSystem) + NbEffParties + year",
  "RC2" = "k_2_pct ~ polarization_parties + factor(ElectSystem) + NbEffParties"
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
  gof_omit = "AIC|BIC|Log.|RMSE", 
  output = "k2ols.tex"
)

### Mundlak-Model (Supplementary Materials):
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
  k_2_pct ~ polarization_within +
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
  coef_map = coef_map,
  output = "Mundlak.tex"
)

###############################################################################################
#### Regressions
#### Mikro-Level
optk_df <- readRDS("DATA/optkdf.RDS") # generated through the cluster function above

# join back to original data
cses.all <- cses.df  %>%
  left_join(x = ., y = optk_df, by = c("case_ID", "ID")) %>% 
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
