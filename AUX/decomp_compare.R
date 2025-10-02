head(resRC.df)
head(ic2res.df)
merged.df <- ic2res.df %>% left_join(x = ., y = resRC.df, by = "id") %>%
  dplyr::select(., c("id", "Gini.within", "Gini.between", "BM.within", "Additive.between"))

ggplot(data = merged.df, 
        aes(x = Gini.within, y = BM.within)) +
  geom_point()


diff.df <- merged.df %>% select(., c("id", "Gini.within", "BM.within")) %>%
  mutate(diff = abs(Gini.within - BM.within))


id193 <- france_theil.df %>% filter(., id == 193)
id1192 <- france_theil.df %>% filter(., id == 1192)
