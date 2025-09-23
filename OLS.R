setwd("~/Documents/Research/Dichotomous/github/Dichotomous/")
library(dplyr)
library(texreg)
library(tidyverse)
library(modelsummary)
library(estimatr)

##
## Grenoble data
idsradrightgre <- read.csv("DATA/idsradrightgre.csv", 
                           sep="", header = T) 
idsradleftgre <- read.csv("DATA/idsradleftgre.csv", 
                           sep="", header = T) 
greno.avg <- read.csv("DATA/OLSgrenoble.csv", header = T) %>%
  mutate(ext.left =  ifelse(id %in% idsradleftgre$x, 1, 0),
         ext.right = ifelse(id %in% idsradrightgre$x, 1, 0)) %>%
  mutate(Age = as.numeric(AGE))
rm(idsradleftgre, idsradrightgre)

ols01 <- lm(formula = phi.value ~ optk2 + Age + 
              GENDER + Educ.lvl + 
              factor(ext.left) + factor(ext.right),
            data = greno.avg)

summary(ols01)
modelsummary(ols01, statistic = "[{conf.low}, {conf.high}]")
#gtsummary::tbl_regression(ols01)
#texreg::texreg(ols01, single.row = T, label = "tb.ols", booktabs = T,
#               custom.model.names = "Grenoble")
#texreg::screenreg(ols01)
#########################################################################
## Graz data
extr.graz <- read.csv("DATA/extrGraz.csv", header = T)
graz.avg <- read.csv("DATA/OLSgraz.csv", header = T) %>%
  left_join(x = ., y = extr.graz, by = "id")
rm(extr.graz)
ols02 <- lm(formula = phi.value ~ optk2 + Age.num + Gender + Educ.lvl +
              factor(ext.left) + factor(ext.right),
            data = graz.avg)
summary(ols02)
modelsummary(ols02, statistic = "[{conf.low}, {conf.high}]")
#gtsummary::tbl_regression(ols02)

# France22 data
extr.fr22 <- read.csv("DATA/extrFR22.csv", header = TRUE)
fra.avg <- read.csv("DATA/OLSfra.csv", header = T) %>% distinct(.) %>%
  left_join(x = ., y = extr.fr22, by = "id")
rm(extr.fr22)
ols03 <- lm(formula = phi.value ~ optk2 +  Age + Gender + Educ.lvl +
              factor(ext.right) + factor(ext.left), 
            data = fra.avg)
summary(ols03)

modelsummary(list(ols01, ols02, ols03),
            stars = TRUE,
            statistic = "[{conf.low}, {conf.high}]", 
            coef_map = "optk2k=2+"
          )



# Common dataset for FE-Regression
grenoble.part <- greno.avg %>%
  mutate(Gender = trimws(GENDER)) %>%
  mutate(Age.num = as.numeric(AGE)) %>%
  select(., c("id", "phi.value", "optk2", "Age.num", "Gender", "Educ.lvl", "ext.right", "ext.left"))

graz.part <- graz.avg %>%
  mutate(Gender = ifelse(Gender == "female", "F", 
                         ifelse(Gender == "male", "M", NA))) %>%
  mutate(Age.num = as.numeric(Age.num)) %>%
  select(., c("id", "phi.value", "optk2", "Age.num", "Gender", "Educ.lvl", "ext.right", "ext.left"))

fra.part <- fra.avg %>%
  mutate(Age.num = as.numeric(Age)) %>%
  mutate(Gender = ifelse(Gender == "h", "M", "F")) %>%
  select(., c("id", "phi.value", "optk2", "Age.num", "Gender", "Educ.lvl", "ext.right", "ext.left"))

ols.df <- rbind(grenoble.part %>% mutate(df = "Grenoble"),
                graz.part %>%     mutate(df = "Graz"), 
                fra.part  %>%     mutate(df = "France")) %>% 
  select(., -c("id")) %>%
  mutate(optk2 = as.factor(optk2),
         Gender = as.factor(Gender),
         df = as.factor(df))

mod01 <- phi.value ~ optk2 +  Age.num + Gender + Educ.lvl
mod02 <- phi.value ~ optk2 +  Age.num + Gender + Educ.lvl + factor(ext.right) + factor(ext.left)

ols.list <- lapply(unique(ols.df$df), FUN = function(i){
  lm(formula = mod02, 
     data = ols.df %>% filter(., df == i) ) 
  } )
texreg::screenreg(l = ols.list)
modelsummary(ols.list, stars = TRUE)


custom.names <- c("(Intercept)" = "Intercept",
                  "optk2B" = "$\tilde{k}>2$",
                  "Age.num" = "Age",
                  "GenderMale" = "Gender (Male)",
                  "Educ.lvl" = "Education Level",
                  "factor(ext.right)1" = "Dummy extr. right",
                  "factor(ext.left)1" =  "Dummy extr. left")  

texreg::texreg(l = ols.list, single.row = T, 
               custom.model.names = c("Grenoble", "Graz", "France 22"),
               booktabs = T, use.package = F, 
#               file = "~/Documents/Research/Dichotomous/tbols.tex",
               custom.coef.names =  custom.names)


## In a robustness check, we merged the three datasets to a common one and ran a
## regression with the survey fixed effects variables. The coefficient is still
# ## negative in this model specification $(-0.08,  [-0.10; -0.06]$) on a 99\% 
# ## significance level, whereas no significant effect is present in any of the controls. 
library(plm)
model_fe <- plm(mod01, 
                data = ols.df, 
                index = "df", 
                model = "within")  # Fixed effects model
summary(model_fe)
lmtest::bptest(model_fe, studentize = TRUE) ## This is important: heteroskadicity is present!
lmtest::coeftest(model_fe, vcov = vcovHC(model_fe, type = "HC1"))  # HC1 
modelsummary(model_fe, 
 # statistic = "{p.value} [{conf.low}, {conf.high}]", 
  vcov = "HC1")   # here we consider the robust SE in the data presentation


model_fe2 <- plm(mod02, 
                data = ols.df, 
                index = "df", 
                model = "within")  # Fixed effects model
summary(model_fe2)
modelsummary(model_fe2)


#gtsummary::tbl_regression(model_fe2, conf.level = 0.95)
#texreg(model_fe2, , single.row = T, 
#       stars = c(0.01),
#                  booktabs = T, use.package = F, 
#                  #               file = "~/Documents/Research/Dichotomous/tbols.tex",
#                  custom.coef.names =  custom.names[-1]
#                  )#


###########################################################
## Section 7 (Soziodemografics)
## Graz
rawdata.graz <- read.csv("regaustria.csv", header = T)
educ.df <- data.frame(Educ = unique(rawdata.graz$Educ),
                      Educ.lvl = c(6, 5, 2, 3, 4, 1, NA))
age.df <- data.frame(Age = sort(unique(rawdata.graz$Age)),
                     Age.num = 1:7)
gender.df <- data.frame(Gender = sort(unique(rawdata.graz$Gender)),
                        Gender.bin = c(NA, "F", "M"))
extr.graz <- read.csv("extrGraz.csv", header = T)
regdata.graz <- rawdata.graz %>% 
  left_join(x = ., y = educ.df, by = "Educ") %>%
  left_join(x = ., y = age.df , by = "Age") %>%
  left_join(x = ., y = gender.df, by = "Gender") %>%
  left_join(x = ., y = extr.graz, by = "id") %>%
  mutate(Age = Age.num, Gender = Gender.bin) %>%
  mutate(k2all = ifelse(id %in% modera.ids.graz, 1, 0),  # from RFGraz.R, line 264f
         k2strong = ifelse(id %in% strong.ids.graz, 1, 0)) %>%
  select(., c("id", "Gender", "Age", "Educ.lvl", 
              starts_with("WDP."), starts_with("k2"), "ext.right", "ext.left")) 
  
rm(rawdata.graz, educ.df, age.df, gender.df, extr.graz)

# Grenoble
rawdata.grenoble <- read.csv("reggrenoble.csv", header = T) %>% 
  mutate(Educ.lvl = as.numeric(ifelse(EDUC  == " S", 3, EDUC)))
age.df <- data.frame(AGE = sort(unique(rawdata.grenoble$AGE)),
                     Age.num = 1:6)
idsradrightgre <- read.csv("~/Documents/Research/Dichotomous/idsradrightgre.csv", 
                           sep="", header = T) 
idsradleftgre <- read.csv("~/Documents/Research/Dichotomous/idsradleftgre.csv", 
                          sep="", header = T) 
regdata.grenoble <- rawdata.grenoble %>% 
  left_join(x = ., y = age.df , by = "AGE") %>%
  select(., c("id", "GENDER", "Age.num", "Educ.lvl", starts_with("WDP."))) %>%
  mutate(Age = Age.num, Gender = GENDER) %>%
  mutate(k2all = ifelse(id %in% modera.ids.grenoble, 1, 0),  # from RFGraz.R, line 264f
         k2strong = ifelse(id %in% strong.ids.grenoble, 1, 0)) %>%
  mutate(ext.left =  ifelse(id %in% idsradleftgre$x, 1, 0),
         ext.right = ifelse(id %in% idsradrightgre$x, 1, 0)) %>%
  select(., colnames(regdata.graz))
rm(rawdata.grenoble, age.df, idsradleftgre, idsradrightgre)

## France
rawdata.france <- read.csv("regfrance.csv", header = T) %>%  
  mutate(studies = na_if(studies, "nspp")) %>%
  mutate(studies = na_if(studies, "")) %>%
  mutate(Educ.lvl = ifelse(studies == "primaire", 1, ifelse(studies == "secondaire", 2, 3))) %>%
  mutate(Gender = na_if(Gender, "nspp")) %>%
  mutate(Gender = na_if(Gender, "")) %>%
  mutate(Age = na_if(Age, "nspp")) %>%
  mutate(Age = na_if(Age, "")) 
age.df <- data.frame(Age = sort(unique(rawdata.france$Age)),
                     Age.num = c(7, 1:6))
extr.fr22 <- read.csv("extrFR22.csv", header = TRUE)
regdata.france <- rawdata.france %>%
  left_join(x = ., 
            y = age.df, by = "Age") %>%
  left_join(x = ., 
            y = extr.fr22, by = "id") %>%
  mutate(Gender = ifelse(Gender == "f", "F", "M")) %>%
  mutate(Age = Age.num) %>%
  mutate(k2all = ifelse(id %in% modera.ids.france, 1, 0),  # from RFGraz.R, line 264f
         k2strong = ifelse(id %in% strong.ids.france, 1, 0)) %>%
  select(., colnames(regdata.graz)) 
rm(rawdata.france, age.df, extr.fr22)

### Regression models
regfun <- function(m){
  logreg <- glm(formula = get(m) ~ Age + Gender + Educ.lvl,
      data = regdata.france,
      family = "binomial")
  return(logreg)
}
m.list <- as.list(colnames(regdata.france[5:9]))
reglist.france <- lapply(m.list, regfun)


## Regressions with Fixed Effects
FE.df <- rbind(
  regdata.grenoble %>% mutate(df = "Grenoble"),
  regdata.graz     %>% mutate(df = "Graz"),
  regdata.france   %>% mutate(df = "France22")
) %>% mutate(Gender = trimws(Gender))

regfun.FE <- function(m){
  logreg <- glm(formula = get(m) ~ Age + Gender + Educ.lvl + factor(df),
                data = FE.df,
                family = "binomial")
  return(logreg)
}
m.list <- as.list(colnames(FE.df[5:9]))
reglist.FE <- lapply(m.list, regfun.FE)

screenreg(reglist.FE, custom.coef.names =  custom.names)

## with interaction effect
FE.sc.df <- FE.df %>%
  mutate(Age.sc = scale(Age, center = T, scale = F),
         Educ.lvl.sc = scale(Educ.lvl, center = T, scale = F))

regfunI.FE <- function(m){
  logreg <- glm(formula = get(m) ~ Age.sc * Educ.lvl.sc + Gender + factor(df),
                data = FE.sc.df,
                family = "binomial")
  return(logreg)
}
m.list <- as.list(colnames(FE.sc.df[5:9]))
reglist.FEI <- lapply(m.list, regfunI.FE)

# with indvidual (ext.right / ext.left) dummies 
custom.names2 <- c("(Intercept)" = "Intercept",
                  "Age.sc" = "Age",
                  "GenderM" = "Gender (Male)",
                  "Educ.lvl.sc" = "Education Level",
                  "factor(ext.right)1" = "Dummy extr. right",
                  "factor(ext.left)1" =  "Dummy extr. left",
                  "factor(df)Grenoble" = "FE-Grenoble",
                  "factor(df)Graz" = "FE-Graz")
regfunD.FE <- function(m){
  logreg <- glm(formula = get(m) ~ Age.sc + Educ.lvl.sc + Gender + 
                  factor(df) + factor(ext.right) + factor(ext.left),
                data = FE.sc.df,
                family = "binomial")
  return(logreg)
}
m.list <- as.list(colnames(FE.sc.df[5:9]))
reglistD.FE <- lapply(m.list, regfunD.FE)
screenreg(l = reglistD.FE)
plotreg(l = reglistD.FE, 
        custom.coef.names = custom.names2)

## with individual (ext.right / ext.left) dummies and interaction:
regfunDI.FE <- function(m){
  logreg <- glm(formula = get(m) ~ Age.sc * Educ.lvl.sc + Gender + 
                  factor(df) + factor(ext.right) + factor(ext.left),
                data = FE.sc.df,
                family = "binomial")
  return(logreg)
}
m.list <- as.list(colnames(FE.sc.df[5:9]))
reglist.FEDI <- lapply(m.list, regfunDI.FE)
screenreg(l = reglist.FEDI)
plotreg(l = reglist.FEDI)

#screenreg(reglist.FEI, ci.force = T, ci.test = 1)
convert_to_or <- function(model) {
  tr <- texreg::extract(model)  # Extract model summary
  tr@coef <- exp(tr@coef)  # Exponentiate coefficients
  tr@ci.low <- exp(tr@ci.low)  # Exponentiate lower CI
  tr@ci.up <- exp(tr@ci.up)  # Exponentiate upper CI
  return(tr)
}

# Apply transformation to all models in reglist.FEI
reglist_or <- lapply(reglist.FEDI, convert_to_or)

# Create LaTeX table with transformed OR values
texreg(l = reglist_or,
       custom.model.names = colnames(FE.df[5:9]),
       caption = "Dichotomy and sociodemografic factors - Interaction term model",
       label = "tb.logregFEI",
 #      file = "InteractlogregFE.tex", 
       ci.test = 1,
       booktabs = T, 
       caption.above = T,
       single.row = F, 
       ci.force = TRUE)

plotreg(l = reglist_or,
        omit.coef = "(Intercept)",
#        type = "forest", 
        ci.test = 1,
        custom.coef.names = custom.names2
) + 
  theme_gray(base_size = 22)


library(car)
vif(reglist.FEI[[1]])
lapply(1:5, function(i){vif(reglist.FEI[[i]])}  )
#library(margins)
#mfx_results <- margins(reglist.FEI[[1]], 
#                       variables = c("Age.sc", "Age.sc:Educ.lvl.sc"),
#                       at = list("Educ.lvl.sc" = c(-2.45, -1.45, -0.45, 0.55, 1.55, 2.55)))
#summary(mfx_results)
#mfx_results <- marginal_effects(model = reglist.FEI[[1]])
#fei.unscaled <- glm(formula = WDP.Theil ~ Age * Educ.lvl + Gender + factor(df),
#                    data = FE.df,
#                    family = "binomial")
#summary(fei.unscaled)
#mfx_results <- marginal_effects(model = fei.unscaled)
#summary(mfx_results)
ggplot(data = FE.df,
       aes(x = Age, colour = df, group = df)) +
  geom_histogram() +
  facet_wrap(~df)

ggplot(data = FE.df,
       aes(x = Age, y = factor(Educ.lvl), group = df, colour = df)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~df)

