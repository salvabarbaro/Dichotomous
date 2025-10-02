## how often is a party not rated 
## Graz
aut.rat <- austria.df %>% dplyr::select(., starts_with("Rat."))

nas <- aut.rat %>%  reframe(across(starts_with("Rat."),
                   ~ mean(is.na(.)),
                   .names = "share_na_{.col}"))

cor.df <- data.frame(party = c("SPÖ", "ÖVP", "FPÖ", "Green", "KPÖ", "NEOS"),
                     match = c(.606, .7, .965, .798, .670, .578),
                     nas = as.numeric(nas) 
)



ggplot(data = cor.df, aes(x = match, y = nas)) + 
  geom_point() + geom_line() +
  theme_gray(base_size = 22) + 
  labs(y = "Share of missing values in the ratings")


# Grenoble
gre.rat <- grenoble.df %>% dplyr::select(., starts_with("EV_"))
nas <- gre.rat %>%  reframe(across(starts_with("EV_"),
                   ~ mean(is.na(.)),
                   .names = "share_na_{.col}"))


cor.df <- data.frame(party = c("NDA", "MLP", "EM", "BH", "NA", "PP", "JC", "JL", "JLM", "FA", "FF"),
                     match = c(.819,  .971,  .838, .832, .669, .704, .764, .713, .843, .809, .932),
                     nas = as.numeric(nas) 
)




ggplot(data = cor.df, aes(x = match, y = nas, label = party)) + 
  geom_point(size = 3) + 
  geom_line() +
  geom_text(vjust = -0.5, hjust = 0.5, size = 6, col = "red") +   # adjust placement and size
  theme_gray(base_size = 22) + 
  labs(y = "Share of missing values in the ratings",
       x = "Matching value")
ggsave("~/Documents/Research/Dichotomous/git/67b5f34c104b85acf4a11317/Reports/popular.pdf", width = 16, height = 9)



