#########################################################################
### Supplementary analyses to Section 4 (Inequality-Measure Approach)
### We analyse how the share of individuals drops depending on the number of approved candidates
# 1. France22
nbapprov.df <- france_long.df %>% group_by(id) %>%
  reframe(id = id, 
          nb.approv = sum(Approval, na.rm = T)   
  ) %>% distinct()
head(nbapprov.df)
france_theil.df <- france_theil.df %>% left_join(x = ., y = nbapprov.df, by = "id")
#
df1 <- france_theil.df %>% filter(., nb.approv <3)
df2 <- france_theil.df %>% filter(., nb.approv == 3)
df3 <- france_theil.df %>% filter(., nb.approv == 4)
df4 <- france_theil.df %>% filter(., nb.approv >4)
df5 <- france_theil.df %>% filter(., nb.approv > 3)


ic2res <- lapply(unique(france_theil.df$id), ic2decomp.fun, 
                  data = df4)
ic2res.df <- ic2res %>% do.call(rbind, .) %>%
  mutate(WDP.Theil = ifelse(Theil.within < Theil.between, 1, 0),
         WDP.Gini =  ifelse(Gini.within < Gini.between, 1, 0),
         WDP.Atkinson = ifelse(Atkinson.within < Atkinson.between, 1, 0),
         WDP.SCV = ifelse(SCV.within < SCV.between, 1, 0)
#         WDP.VAR = ifelse(VAR.within < VAR.between, 0, 1)
#         WDP.MLD = ifelse(MLD.within < MLD.between, 1, 0)
        )
compute_wdp_shares(ic2res.df)
#
################################################################################
# 2. Grenoble
nbapprov.df <- grenoble_theil.df %>% group_by(id) %>%
  reframe(id = id, 
          nb.approv = sum(Approval, na.rm = T)   
  ) %>% distinct()

head(nbapprov.df)
table(nbapprov.df$nb.approv)
grenoble_theil.df <- grenoble_theil.df %>% left_join(x = ., y = nbapprov.df, by = "id")
#
df1 <- grenoble_theil.df %>% filter(., nb.approv <3)
df2 <- grenoble_theil.df %>% filter(., nb.approv == 3)
df3 <- grenoble_theil.df %>% filter(., nb.approv == 4)
df4 <- grenoble_theil.df %>% filter(., nb.approv >4)


ic2res <- lapply(unique(grenoble_theil.df$id), ic2decomp.fun, 
                  data = df4)
ic2res.df <- ic2res %>% do.call(rbind, .) %>%
  mutate(WDP.Theil = ifelse(Theil.within < Theil.between, 1, 0),
         WDP.Gini =  ifelse(Gini.within < Gini.between, 1, 0),
         WDP.Atkinson = ifelse(Atkinson.within < Atkinson.between, 1, 0),
         WDP.SCV = ifelse(SCV.within < SCV.between, 1, 0)
#         WDP.VAR = ifelse(VAR.within < VAR.between, 0, 1)
#         WDP.MLD = ifelse(MLD.within < MLD.between, 1, 0)
        )
compute_wdp_shares(ic2res.df)
# 3. Graz
nbapprov.df <- austria_long.df %>% group_by(id) %>%
  reframe(id = id, 
          nb.approv = sum(Approval, na.rm = T)   
  ) %>% distinct()
head(nbapprov.df)
austria_theil.df <- austria_theil.df %>% left_join(x = ., y = nbapprov.df, by = "id")
table(austria_theil.df$nb.approv)
#
df1 <- austria_theil.df %>% filter(., nb.approv <3)
df2 <- austria_theil.df %>% filter(., nb.approv == 3)
df3 <- austria_theil.df %>% filter(., nb.approv == 4)
df4 <- austria_theil.df %>% filter(., nb.approv >4)

ic2res <- lapply(unique(austria_theil.df$id), ic2decomp.fun, 
                  data = df4)
ic2res.df <- ic2res %>% do.call(rbind, .) %>%
  mutate(WDP.Theil = ifelse(Theil.within < Theil.between, 1, 0),
         WDP.Gini =  ifelse(Gini.within < Gini.between, 1, 0),
         WDP.Atkinson = ifelse(Atkinson.within < Atkinson.between, 1, 0),
         WDP.SCV = ifelse(SCV.within < SCV.between, 1, 0)
#         WDP.VAR = ifelse(VAR.within < VAR.between, 0, 1)
#         WDP.MLD = ifelse(MLD.within < MLD.between, 1, 0)
        )
compute_wdp_shares(ic2res.df)

rm(df1, df2, df3, df4)

nbapp.df <- data.frame(
  Grenoble   = c(0.595, .777, .849, .822),
  Graz       = c(.409, .752, .758, .750), 
  France22   = c(.052, .589, .789, .788), 
  France22NA = c(.4476, .730, .859, .848)
)

# Define your x-axis categories
nbapp.df$Index <- factor(c("{1,2}", "3", "4", ">4"),
                         levels = c("{1,2}", "3", "4", ">4"))

# Pivot to long format
nbapp.long <- nbapp.df %>%
  pivot_longer(cols = c(Grenoble, Graz, France22, France22NA),
               names_to = "Dataset",
               values_to = "Value")

# Plot
ggplot(nbapp.long, aes(x = Index, y = Value, colour = Dataset, group = Dataset)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Number of Approved Candidates", y = "Value", colour = "Dataset") +
  theme_gray(base_size = 22) +
  labs(y = "Share", x = "Number of approved alternatives")
ggsave("~/Documents/Research/Dichotomous/git/67b5f34c104b85acf4a11317/nbapproved.pdf", width = 16, height = 9)
