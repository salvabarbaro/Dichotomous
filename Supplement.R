## Supplementary material
## Descriptive Statistics

plot.data <- grenoble.df %>% select(., c("id", starts_with("EV"))) %>%
  pivot_longer(cols = starts_with("EV")) %>%
  mutate(name = gsub("EV_", "", name))


approv.df <- data.frame(Candidates = colnames(grenoble.df %>% select(starts_with("AV"))),
                        approv.count = sapply(grenoble.df %>% select(., starts_with("AV")), sum)) %>%
  mutate(approv.share = approv.count / nrow(grenoble.df)) %>%
  mutate(name = gsub("AV_", "", Candidates)) %>%
  mutate(plurality = c(0.012, 0.024, 0.351, 0.151, 0.003, .005, 0.001, 0.004, .323, .004, .086))



plot.data2 <- plot.data %>% 
  left_join(x = ., y = approv.df, by = "name")
bp.grenoble <- ggplot(data = plot.data2,
       aes(y = value, group = name, colour = name)) +
  geom_boxplot(aes(x = name)) +
  geom_text(aes(x = name, y = 0.5, 
                label = round(approv.share,2)), 
            colour = "black", size = 6) +
  geom_text(aes(x = name, y = 0.25, 
                label = round(plurality, 2)), 
            colour = "red", size = 6) +
  theme_gray(base_size = 22) + theme(legend.position = "none") +
  labs(y = "rating values", x = "Candidate")
ggsave(plot = bp.grenoble, "boxplotGrenoble.pdf", width = 16, height = 9)
##################################################################################
ratings.df <- austria.df %>% select(starts_with("Rat.")) %>%
  pivot_longer(cols = everything(), 
               names_to = "Variable", 
               values_to = "Value") %>%
  mutate(Variable = gsub("Rat.", "", Variable))

bp.graz <- ggplot(ratings.df, aes(x = Variable, y = Value, fill = Variable)) +
  geom_boxplot(notch = FALSE) +
  scale_fill_manual(values = c(
    "FPÖ" = "chocolate",
    "SPÖ" = "red",
    "ÖVP" = "gray8",
    "Green" = "forestgreen",
    "NEOS" = "yellow2",
    "KPÖ" = "purple" )) + 
  labs(title = "Ratings", x = "", y = "Value") +
  theme_gray(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave("boxplotGraz.pdf", width = 16, height = 9, plot = bp.graz)
##
sapply(austria.df %>% select(., starts_with("Approv.")), sum, na.rm = T) / 937

trichfun <- function(p) {
  tv <- table(austria.df[[p]]) / length(austria.df[[p]][!is.na(austria.df[[p]])])
  return(tv)
}

trich.list <- lapply(colnames(austria.df)[17:22], trichfun)
trich.df <- do.call(rbind, trich.list) %>% as.data.frame(.) %>%
  mutate(temp = colnames(austria.df)[17:22], .before = "-1") %>%
  mutate(Party = gsub("Trich.", "", temp), .before = "temp") %>%
  select(., -c("temp"))

stargazer::stargazer(trich.df, summary = F, rownames = F)


###########################################################
## France22
plot.data <- france22.df %>% select(., c("id", starts_with("EV"))) %>%
  pivot_longer(cols = starts_with("EV")) %>%
  mutate(name = gsub("EV_", "", name))


approv.df <- data.frame(Candidates = colnames(france22.df %>% select(starts_with("AV"))),
                        approv.count = sapply(france22.df %>% select(., starts_with("AV")), sum)) %>%
  mutate(approv.share = approv.count / nrow(france22.df)) %>%
  mutate(name = gsub("AV_", "", Candidates)) %>%
  mutate(plurality = c(.001, .008, .008, .023, .095, .022, .022, .012, .017, .121, .661, .010))

#nas <- c("", "abstention", "blanc", "m18", "ni", "nspp", "incertitude")
#plurality <- ifelse(france22.df$plurality %in% nas, NA, france22.df$plurality)
#table(plurality)/length(plurality[is.na(plurality)==F])


plot.data2 <- plot.data %>% 
  left_join(x = ., y = approv.df, by = "name")
bp.france <- ggplot(data = plot.data2,
                      aes(y = value, group = name, colour = name)) +
  geom_violin(aes(x = name)) +
  geom_text(aes(x = name, y = 0.5, 
                label = round(approv.share,2)), 
            colour = "black", size = 6) +
  geom_text(aes(x = name, y = 0.25, 
                label = round(plurality*100, 3)), 
            colour = "red", size = 6) +
  theme_gray(base_size = 22) + theme(legend.position = "none") +
  labs(y = "rating values", x = "Candidate")
ggsave(plot = bp.france, "boxplotFrance.pdf", width = 16, height = 9)




