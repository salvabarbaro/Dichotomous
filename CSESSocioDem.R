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
mod01.sat <- "bin.k2 ~ SatisfactionDem + Age + Gender + Education + IncomeQ | case_ID"
mod01.both <- "bin.k2 ~ Dist0 + SatisfactionDem + Age + Gender + Education + IncomeQ | case_ID"


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
  coef_map = c(
    "IncomeQ" = "Income",
    "GenderM" = "Male",
        "Age" = "Age",
    "SatisfactionDem" = "Satisf. Demo.",
    "Dist0" = "Ideology")
) + theme_bw(base_size = 24) + geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_viridis_d()
ggsave("RegressionPlots.pdf", width = 16, height = 9)


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
