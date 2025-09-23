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
logregmod.WDP <- "WDP.Gini ~ Age + Gender + Educ.lvl + factor(df)"
logregmod.WDP <- "k2strong ~ Age + Gender + Educ.lvl + factor(df)"
logregmod.RC1  <- "WDP.Gini ~ Age + Gender + Educ.lvl + ext.right + ext.left"
logreg.main   <- glm(formula = logregmod.WDP, 
                      family = "binomial", 
                      data = logregdata.df)
modelsummary(logreg.main, exponentiate = T,
              statistic = "[{conf.low}, {conf.high}]")



regfun <- function(m){
  logreg <- glm(formula = get(m) ~ Age + Gender + Educ.lvl,
      data = regdata.fran,
      family = "binomial")
  return(logreg)
}
m.list <- as.list(colnames(regdata.fran)[5:9])
reglist.france <- lapply(m.list, regfun)
modelsummary(reglist.france, exponentiate = T)
