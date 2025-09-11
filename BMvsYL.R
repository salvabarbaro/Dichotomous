b.df <- ic2res.df %>% select(., c("id", "Gini.between")) %>%
  left_join(x = ., 
            y = ic2resBM.df %>% select(., c("id", "Gini.between")),
            by = "id" 
          ) %>%
  dplyr::rename(YL = Gini.between.x,
                BM = Gini.between.y)

ggplot(data = b.df, aes(x = YL, y = BM)) +
    geom_point()


