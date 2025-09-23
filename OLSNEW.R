### OLS NEW  
## Replication Files for the OLS Regresssions
## Generates Table 8 and Fig. 5. 
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

#############################################################################
## I. Data preparation (DP)
#  I.1 DP-Grenoble
idsradrightgre <- read.csv("DATA/idsradrightgre.csv", 
                           sep="", header = T) 
idsradleftgre <- read.csv("DATA/idsradleftgre.csv", 
                           sep="", header = T) 
greno.avg <- read.csv("DATA/OLSgrenoble.csv", header = T) %>%
  mutate(ext.left =  ifelse(id %in% idsradleftgre$x, 1, 0),
         ext.right = ifelse(id %in% idsradrightgre$x, 1, 0)) %>%
  mutate(Age = as.numeric(AGE))
rm(idsradleftgre, idsradrightgre)
# I.2 DP-Graz
extr.graz <- read.csv("DATA/extrGraz.csv", header = T)
graz.avg <- read.csv("DATA/OLSgraz.csv", header = T) %>%
  left_join(x = ., y = extr.graz, by = "id")
rm(extr.graz)
# I.3 DP-France22
extr.fr22 <- read.csv("DATA/extrFR22.csv", header = TRUE)
fra.avg <- read.csv("DATA/OLSfra.csv", header = T) %>% distinct(.) %>%
  left_join(x = ., y = extr.fr22, by = "id")
rm(extr.fr22)
# Common dataset 
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
#################################################################################
## II. Main Model Specification
mod.main <- "phi.value ~ optk2 + Age.num + Gender + Educ.lvl"
Main.Grenoble <- lm(formula = mod.main,
                    data = ols.df %>% dplyr::filter(., df == "Grenoble"))
Main.Graz <-     lm(formula = mod.main,
                    data = ols.df %>% dplyr::filter(., df == "Graz"))
Main.France22 <- lm(formula = mod.main,
                    data = ols.df %>% dplyr::filter(., df == "France"))
# II.Tab (Table 8 in the paper)
modelsummary(list(Main.Grenoble, Main.Graz, Main.France22),
            stars = TRUE,
            statistic = "[{conf.low}, {conf.high}]")
####################################################################################
## II. RC1 - Same as Main, but with left/right dummies
mod.RC1 <- "phi.value ~ optk2 + Age.num + Gender + Educ.lvl + factor(ext.left) + factor(ext.right)"
RC1.Grenoble <- lm(formula = mod.RC1,
                    data = ols.df %>% dplyr::filter(., df == "Grenoble"))
RC1.Graz <-     lm(formula = mod.RC1,
                    data = ols.df %>% dplyr::filter(., df == "Graz"))
RC1.France22 <- lm(formula = mod.RC1,
                    data = ols.df %>% dplyr::filter(., df == "France"))
# II.Tab (Table 8 in the paper)
modelsummary(list(RC1.Grenoble, RC1.Graz, RC1.France22),
            stars = TRUE,
            statistic = "[{conf.low}, {conf.high}]")
######################################################################################
## III. RC2 - Same as Main, with robust SE  # alternatively: run lm_robust()
modelsummary(list(Main.Grenoble, Main.Graz, Main.France22),
            stars = TRUE,
            statistic = "[{conf.low}, {conf.high}]",
            vcov = "HC1")
#######################################################################################
## IV. Common dataset, Fixed effects model, Main
RC3 <-  plm(mod.main, 
            data = ols.df, 
            index = "df", 
            model = "within")  
modelsummary(RC3, stars = TRUE, statistic = "[{conf.low}, {conf.high}]", vcov = "HC1")
########################################################################################
## V. Common dataset, Fixed effects model, With left/right dummies
RC4 <-  plm(mod.RC1, 
            data = ols.df, 
            index = "df", 
            model = "within")  
modelsummary(RC4, stars = TRUE, statistic = "[{conf.low}, {conf.high}]", vcov = "HC1")
##########################################################################################
## VI. Result visualization
## VI.1 Collect the models
mods <- list(
  "Grenoble — Main"            = Main.Grenoble,            # classical SE
  "Grenoble — RC1 (+L/R)"      = RC1.Grenoble,             # classical SE 
  "Grenoble — RC2 (robust)"    = Main.Grenoble,            # HC1
  "Graz — Main"                = Main.Graz,                # classical SE
  "Graz — RC1 (+L/R)"          = RC1.Graz,                 # classical SE
  "Graz — RC2 (robust)"        = Main.Graz,                # HC1
  "France22 — Main"              = Main.France22,            # classical SE
  "France22 — RC1 (+L/R)"        = RC1.France22,             # classical SE
  "France22 — RC2 (robust)"      = Main.France22,            # HC1
  "Common — FE (RC3)"          = RC3,                      # cluster by group
  "Common — FE + L/R (RC4)"    = RC4                       # cluster by group
)
## VI.2 vcov-assignement:
##    - OLS: HC1
##    - FE (plm): HC1 clustered by panel unit (“group”)
vc_list <- mapply(function(nm, m) {
  if (grepl("RC2 \\(robust\\)", nm)) {
    "HC1"  # robust (non-clustered!) for RC2
  } else if (inherits(m, "plm")) {
    function(x) vcovHC(x, type = "HC1", cluster = "group")  # FE: cluster by panel unit
  } else {
    NULL   # classical SE for Main and RC1
  }
}, names(mods), mods, SIMPLIFY = FALSE)

## VI.3 ggplot 
p <- modelplot(
  mods,
  vcov       = vc_list,
  coef_map   = c("optk2k=2+" = "Optimal k=2"),
  conf_level = 0.95
) +
  labs(x = "Estimate (95% CI)", y = NULL) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_minimal(base_size = 22) +
  coord_flip()
p

###################################################################################
## For Fig. 5, we adjusted the plot with customized colours. 
###########################################################
## Extended ggplot (with different colours)
## 1) color map 
cols <- c(
  "Grenoble — Main"            = brewer.pal(9,"Reds")[5],
  "Grenoble — RC1 (+L/R)"      = brewer.pal(9,"Reds")[7],
  "Grenoble — RC2 (robust)"    = brewer.pal(9,"Reds")[9],

  "Graz — Main"                = brewer.pal(9,"Blues")[7],
  "Graz — RC1 (+L/R)"          = brewer.pal(9,"Blues")[8],
  "Graz — RC2 (robust)"        = brewer.pal(9,"Blues")[9],

  "France22 — Main"              = brewer.pal(9,"Greens")[7],
  "France22 — RC1 (+L/R)"        = brewer.pal(9,"Greens")[8],
  "France22 — RC2 (robust)"      = brewer.pal(9,"Greens")[9],

  "Common — FE (RC3)"          = brewer.pal(9,"Purples")[8],
  "Common — FE + L/R (RC4)"    = brewer.pal(9,"Oranges")[7]
)

## 2) Get the tidy data from modelplot (so we fully control aesthetics)
td <- modelplot(
  mods,
  vcov       = vc_list,
  coef_map   = c("optk2k=2+" = "Optimal k=2"),  
  conf_level = 0.95,
  draw       = FALSE
)

## Keep ordering identical to your `mods` list
td$model <- factor(td$model, levels = names(mods))

## 3) Build the plot with manual colours & flipped coords
p2 <- ggplot(td, aes(y = model, x = estimate,
                    xmin = conf.low, xmax = conf.high,
                    colour = model)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(height = 0.18, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = "Estimate (95% CI)", y = NULL) +
  theme_gray(base_size = 22) +
  theme(
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(), legend.position = "bottom") +
  coord_flip()

p2

ggsave("~/Documents/Research/Dichotomous/git/67b5f34c104b85acf4a11317/OLSRC.pdf", width = 16, height = 9, plot = p2)



