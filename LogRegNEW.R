## Replication files for the logistic regression in the Sociodemographics-Section
setwd("~/Documents/Research/Dichotomous/github/Dichotomous/")
library(dplyr)
#library(texreg)
library(tidyverse)
library(modelsummary)
#library(estimatr) ## if lm_robust 
library(ggplot2)
library(sandwich)
library(plm)
library(RColorBrewer)

##################################################################################
## I. Get the data
## Hint: the RF*.R files generate the information regarding affinity to radical forces and the IDs
## associated with WDP, k2_strong, k2_weak.
regdata.gren <- read.csv("DATA/regdataGrenoble.csv", header = T)
regdata.graz <- read.csv("DATA/regdataGraz.csv", header = T)
regdata.fran <- read.csv("DATA/regdataFrance22.csv", header = T)
logregdata.df <- rbind(regdata.gren %>% mutate(df = "Grenoble"),
                       regdata.graz %>% mutate(df = "Graz"), 
                       regdata.fran %>% mutate(df = "France22"))
#######################################################################################
# models
logregmod.WDP <- "WDP.Gini ~ Age + Gender + Educ.lvl"
logregmod.k2 <- "k2strong ~ Age + Gender + Educ.lvl"
logregmod.RC1  <- "WDP.Gini ~ Age + Gender + Educ.lvl + ext.right + ext.left"
logregmod.k2R <- "k2strong ~ Age + Gender + Educ.lvl + ext.right + ext.left"
# Regressions
logreg.WDP.Gren <- glm(formula = logregmod.WDP, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Grenoble"))
#
logreg.WDP.Gra <- glm(formula = logregmod.WDP, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Graz"))
#
logreg.WDP.Fra <- glm(formula = logregmod.WDP, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "France22"))
#
logreg.RC1.Gren <- glm(formula = logregmod.RC1, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Grenoble"))
#
logreg.RC1.Gra <- glm(formula = logregmod.RC1, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Graz"))
#
logreg.RC1.Fra <- glm(formula = logregmod.RC1, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "France22"))
#
logreg.k2.Gren <- glm(formula = logregmod.k2, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Grenoble"))
#
logreg.k2.Gra <- glm(formula = logregmod.k2, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Graz"))
#
logreg.k2.Fra <- glm(formula = logregmod.k2, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "France22"))
#
logreg.k2R.Gren <- glm(formula = logregmod.k2R, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Grenoble"))
#
logreg.k2R.Gra <- glm(formula = logregmod.k2R, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "Graz"))
#
logreg.k2R.Fra <- glm(formula = logregmod.k2R, 
                       family = "binomial", 
                       data = logregdata.df %>% filter(., df == "France22"))
##################################################################################


## Results
mods <- list(logreg.WDP.Gren, logreg.WDP.Gra, logreg.WDP.Fra,
             logreg.RC1.Gren, logreg.RC1.Gra, logreg.RC1.Fra,
             logreg.k2.Gren, logreg.k2.Gra, logreg.k2.Fra,
             logreg.k2R.Gren, logreg.k2R.Gra, logreg.k2R.Fra)

mods <- list("WDP.Main - Grenoble" =logreg.WDP.Gren, 
             "WDP.Main - Graz" = logreg.WDP.Gra, 
             "WDP.Main - France22" = logreg.WDP.Fra,
             "WDP.Extr - Grenoble" = logreg.RC1.Gren, 
             "WDP.Extr - Graz" = logreg.RC1.Gra, 
             "WDP.Extr - France22" =logreg.RC1.Fra,
             "K2.Main - Grenoble" = logreg.k2.Gren, 
             "K2.Main - Graz" =logreg.k2.Gra, 
             "K2.Main - France22" =logreg.k2.Fra,
             "K2.Extr - Grenoble" = logreg.k2R.Gren, 
             "K2.Extr - Graz" =logreg.k2R.Gra, 
             "K2.Extr - France22" =logreg.k2R.Fra)


modelsummary(mods, 
            exponentiate = T,
             statistic = "[{conf.low}, {conf.high}]",
            stars = T)

modelplot(mods[1:6], 
          exponentiate = T, 
          coef_map   = c("Educ.lvl" = "Educ. lvl."),  
          conf_level = 0.9) + 
  geom_vline(xintercept = 1, col = "black") +
  coord_flip() + theme_gray(base_size = 22)




## Nächste Schritte: Neue Variable für Anzahl bewerteter und Anzahl approved candidate.
## K2 auf den CSES-Datensatz anwenden. 












regfun <- function(m){
  logreg <- glm(formula = get(m) ~ Age + Gender + Educ.lvl,
      data = regdata.fran,
      family = "binomial")
  return(logreg)
}
m.list <- as.list(colnames(regdata.fran)[5:9])
reglist.france <- lapply(m.list, regfun)
modelsummary(reglist.france, exponentiate = T)
