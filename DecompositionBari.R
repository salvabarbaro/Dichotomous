library(ineq)
library(dineq)
#library(haven)
library(dplyr)
library(tidyr)

#Functions

{

gini_coefficient <- function(x) {
    
    return(Gini(x, corr = FALSE, na.rm = TRUE))
  }
  
gini_within_add_path <- function(incomes, group) {
    
    total_income <- sum(incomes)
    total_n <- length(incomes)
    
    # Initialize the weighted sum of Gini coefficients
    weighted_gini_add <- 0
    weighted_gini_path <- 0
    
    # Get the unique groups
    unique_groups <- unique(group)
    
    # For each group, calculate its Gini coefficient and weight, and update the weighted sum
    for (g in unique_groups) {
      # Subset incomes and group for the current group
      group_incomes <- incomes[group == g]
      
      # Calculate the Gini coefficient for the current group
      gini_g <- gini_coefficient(group_incomes)
      
      # Calculate the number of elements in the group and the group mean
      n_g <- length(group_incomes)
      q_g <- sum(group_incomes)
      
      # Calculate the weight for the group
      weight_g_add <- (n_g / total_n) * (q_g / total_income)
      
      weight_g_path <- (n_g / total_n)^2
      
      
      # Add the weighted Gini coefficient to the sum
      weighted_gini_add <- weighted_gini_add + weight_g_add * gini_g
      weighted_gini_path <- weight_g_path + weight_g_path * gini_g
      
    }
    
    values <- c(weighted_gini_add, weighted_gini_path)
    # Return the weighted sum of Gini coefficients
    return(values)
  }
  
gini_between_add <- function(incomes, group){
    # Create a data frame with incomes and their respective group
    data <- data.frame(incomes = incomes, group = group)
    
    # Calculate the mean income for each group
    group_means <- tapply(data$incomes, data$group, mean)
    
    # Replace each income with the mean of the group it belongs to
    data$incomes <- group_means[data$group]
    
    #Retuns its Gini coefficient
    return(gini_coefficient(data$incomes))
  } 
  
  
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




#Decompositions over time only Ethnicity

{
  
    
    # Read the dataset
    data <- read.csv(file_name)
    
    
    # Replace values in the Ethnicity column
    data$Ethnicity[data$Ethnicity != 1] <- 2
    
    
    G <-  gini_coefficient(data$income)
    
    values <- gini_within_add_path(data$income, data$Ethnicity)
    
    GwA <- values[1]
    GwP <- values[2]
    
    GbA <- gini_between_add(data$income, data$Ethnicity)
    
    R <- G - GwA - GbA
    
    GbP <- G - GwP - R
    
    R_star <- G - GwP - GbA
    
     
    # Store the results in a dataframe
    results_df <- data.frame(
      Year = year,
      Gini = G,
      G_bet_add = GbA,
      G_wit_add = GwA,
      G_bet_path = GbP,
      G_wit_path = GwP,
      Res_add = R,
      Res_new = R_star
    )
    
   
    
  }
  


set.seed(55234)
df <- data.frame(inc = sample(1000:5000, 100, replace = T),
      group = sample(1:2, 100, replace = TRUE) )
#
# 1. Gini coefficient
gini_coefficient(df$inc)

gini_within_add_path(incomes = df$inc, group = df$group)

## france_theil.df

test.df <- france_theil.df %>% filter(., id == 6)
gini_within_add_path(incomes = test.df$Rating, group = test.df$Approval)[1]
gini_between_add(incomes = test.df$Rating, group = test.df$Approval)



## adjusted script
{
# Read the dataset
    data <- france_theil.df %>% filter(., id == 2215)
    G <-  gini_coefficient(data$Rating)
    values <- gini_within_add_path(data$Rating, data$Approval)
    GwA <- values[1]
    GwP <- values[2]
    GbA <- gini_between_add(data$Rating, data$Approval)
    R <- G - GwA - GbA
    GbP <- G - GwP - R
    R_star <- G - GwP - GbA
    # Store the results in a dataframe
    results_df <- data.frame(
#      Year = year,
      Gini = G,
      G_bet_add = GbA,
      G_wit_add = GwA,
      G_bet_path = GbP,
      G_wit_path = GwP,
      Res_add = R,
      Res_new = R_star
    )
  }
###################
results_df



df <- tibble::tibble(
  Approval = c(0,0,0,0,1,0,1,0,0,0,1,0),
  Rating   = c(2.5,2,2,2.5,3,2,3,2.5,2,2.5,3,2.5)
)

gini_between_add(incomes = df$Rating, group = df$Approval)







#### What we do: we take the Gini.between and compare it with the G_wit_path.
# First step, get a list of G_with_path

G_wit_path_only <- function(y, g) gini_within_add_path(y, g)[2]

res_gwit_path <- france_theil.df %>%
  filter(is.finite(Rating), !is.na(Approval)) %>%   # clean inputs
  group_by(id) %>%
  summarise(G_wit_path = G_wit_path_only(Rating, Approval), .groups = "drop")

res_gwit_path

## Fleurbaey et.al. 
fleurbaey <- ic2res.df %>% dplyr::select(., c("id", "Gini.between")) %>%
  left_join(x = ., y = res_gwit_path, by = "id") %>%
  mutate(WDP.Fleurbaey = ifelse(G_wit_path < Gini.between, 1, 0))



G_wit_path_only <- function(y, g) gini_within_add_path(y, g)[2]

res_ids <- france_theil.df %>%
  filter(is.finite(Rating), !is.na(Approval)) %>%
  group_by(id) %>%
  summarise(
    G_wit_path = G_wit_path_only(Rating, Approval),
    G_bet_add  = gini_between_add(Rating, Approval),
    .groups = "drop"
  )

#######
# Pure vector-in / scalar-out BM between component:
bm_between_vec <- function(y, g) {
  y <- as.numeric(y); g <- as.factor(g)
  N  <- length(y); if (N < 2) return(0)
  mu <- mean(y)
  mu_g <- tapply(y, g, mean)
  n_g  <- tapply(y, g, length)
  p    <- n_g / N
  # (1/(2*mu)) * sum_k sum_h p_k p_h |mu_k - mu_h|
  sum(outer(p, p) * abs(outer(mu_g, mu_g, `-`))) / (2 * mu)
}

res_ids <- france_theil.df %>%
  filter(is.finite(Rating), !is.na(Approval)) %>%
  group_by(id) %>%
  group_modify(~{
    y <- .x$Rating
    g <- .x$Approval
    GwP <- gini_within_add_path(y, g)[2]          # within (path)
    GbA <- bm_between_vec(y, g)                   # between (add), safe
    G   <- gini_coefficient(y)
    tibble(
      gini_wit_path = GwP,
      Res_new       = G - GwP - GbA               # R* = G - GwP - GbA
    )
  }) %>%
  ungroup()

#res_ids
fleurbaey <- ic2res.df %>% dplyr::select(., c("id", "Gini.between")) %>%
  left_join(x = ., y = res_ids, by = "id") %>%
  mutate(WDP.Fleurbaey = ifelse(gini_wit_path < (Gini.between - Res_new), 1, 0))
table(fleurbaey$WDP.Fleurbaey)
summary(fleurbaey$Res_new)



###########################################################################################
### Fleurbay-Function

# --- helper: pure vector-in / scalar-out BM between component ---
bm_between_vec <- function(y, g) {
  y <- as.numeric(y); g <- as.factor(g)
  N  <- length(y); if (N < 2) return(0)
  mu <- mean(y)
  mu_g <- tapply(y, g, mean)
  n_g  <- tapply(y, g, length)
  p    <- n_g / N
  # (1/(2*mu)) * sum_k sum_h p_k p_h |mu_k - mu_h|
  sum(outer(p, p) * abs(outer(mu_g, mu_g, `-`))) / (2 * mu)
}

# --- main wrapper ---
# data: your data frame
# id_col, income_col, group_col: column names (strings)
# join_df: optional external df with per-id Gini.between
# join_id_col, join_between_col: column names in join_df
fleurbaey_pipeline <- function(
  data,
  id_col        = "id",
  income_col    = "Rating",
  group_col     = "Approval",
  join_df       = NULL,
  join_id_col   = "id",
  join_between_col = "Gini.between"
) {
  # small safe-call helpers (return NA on error)
  safe_gwP  <- function(y, g) tryCatch(gini_within_add_path(y, g)[2], error = function(e) NA_real_)
  safe_gini <- function(y)    tryCatch(gini_coefficient(y),         error = function(e) NA_real_)

  # compute per-id: gini_wit_path and Res_new = G - GwP - GbA
  res_ids <- data %>%
    filter(is.finite(.data[[income_col]]), !is.na(.data[[group_col]])) %>%
    group_by(.data[[id_col]]) %>%
    group_modify(~{
      y <- .x[[income_col]]
      g <- .x[[group_col]]
      GwP <- safe_gwP(y, g)               # within (path)
      GbA <- bm_between_vec(y, g)         # between (add), safe
      G   <- safe_gini(y)                 # overall Gini
      tibble(
        gini_wit_path = GwP,
        Res_new       = G - GwP - GbA     # R* = G - GwP - GbA
      )
    }) %>%
    ungroup() %>%
    rename(id = !!id_col)

  # If no external join requested, return early
  if (is.null(join_df)) {
    return(list(res_ids = res_ids))
  }

  # Join with external per-id Gini.between and compute WDP.Fleurbaey
  fleurbaey <- join_df %>%
    select(!!join_id_col, !!join_between_col) %>%
    rename(id = !!join_id_col, Gini_between = !!join_between_col) %>%
    left_join(res_ids, by = "id") %>%
    mutate(WDP.Fleurbaey = ifelse(gini_wit_path < (Gini_between - Res_new), 1L, 0L))

  # Summary outputs (as requested in your snippet)
  wdp_table   <- table(fleurbaey$WDP.Fleurbaey, useNA = "ifany")
  resnew_summ <- summary(fleurbaey$Res_new)

  list(
    res_ids      = res_ids,      # two columns: id, gini_wit_path & Res_new
    fleurbaey    = fleurbaey,    # joined df with WDP.Fleurbaey
    wdp_table    = wdp_table,    # table(WDP.Fleurbaey)
    resnew_summ  = resnew_summ   # summary(Res_new)
  )
}
## apply the functions
# Just the per-id results (no join):
out1 <- fleurbaey_pipeline(france_theil.df)
out1$res_ids

# With your external table ic2res.df (must have columns id and Gini.between):
out2 <- fleurbaey_pipeline(
  data = france_theil.df,
  join_df = ic2res.df,
  id_col = "id",
  income_col = "Rating",
  group_col = "Approval",
  join_id_col = "id",
  join_between_col = "Gini.between"
)

out2$res_ids
out2$fleurbaey
out2$wdp_table
out2$resnew_summ
