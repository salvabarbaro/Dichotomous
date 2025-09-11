################################################################################################
############# AUX-Functions ####################################################################
# This script contains functions required to run the Weak-Dichotomous Preference Analyses
# 1. cv_decomp_bm: Decomposes the squared coeff of var with the Bhattacharya–Mahalanobis (BM)-approach
# 2. gini_within_add_path is the within part of the Foster and Shneyerov approach
# 3. gini_coefficient is the simple Gini calculation and equivalent to ineq::ineq(x, type = "Gini")
# 4. bm_between_vec calculates the Gini-between dispersion based on BM. Equivalent to IC2::decompSGini(... decomp = "BM")
# 5. fleurbaey_pipeline calculates the comparision proposed by Fleurbaey et. al. 2025, SocChWelf
# 6. ic2decomp.fun calculates the four different decompositions shown in Fig 3 of the paper
########################################################################################################
# 1. 
# Decompose by squared coefficient of variation (SCV = Var / mu^2)
# using the Bhattacharya–Mahalanobis variance decomposition.
# income = column name with incomes (e.g., "Rating")
# group  = column name with groups   (e.g., "Approval")
cv_decomp_bm <- function(df, income, group, na.rm = TRUE, per_group = FALSE) {
  # pull columns
  y <- df[[income]]
  g <- df[[group]]
  
  # handle NAs
  if (na.rm) {
    keep <- is.finite(y) & !is.na(g)
    y <- y[keep]; g <- g[keep]
  }
  
  N  <- length(y)
  if (N < 2) stop("Need at least two observations.")
  mu <- mean(y)
  
  # totals
  total_var <- var(y)
  total_scv <- total_var / (mu^2)
  
  # group stats
  mu_g  <- tapply(y, g, mean)
  n_g   <- tapply(y, g, length)
  var_g <- tapply(y, g, var)
  var_g[is.na(var_g)] <- 0  # variance is 0 if group size = 1
  p_g   <- n_g / N
  scv_g <- var_g / (mu_g^2)
  
  # BM variance decomposition
  within_var  <- sum(p_g * var_g)
  between_var <- sum(p_g * (mu_g - mu)^2)
  
  # convert to SCV components (exact additivity holds for SCV)
  within_scv  <- within_var  / (mu^2)
  between_scv <- between_var / (mu^2)
  
  # equivalent “weighted group SCV” representation of within SCV
  # sum_g p_g * (mu_g/mu)^2 * scv_g == within_scv
  within_scv_alt <- sum(p_g * (mu_g / mu)^2 * scv_g)
  
  out <- list(
    total_cv        = sqrt(total_scv),   # overall CV
    total_scv       = total_scv,         # exact target being decomposed
    within_scv      = within_scv,
    between_scv     = between_scv,
    check_scv       = within_scv + between_scv,  # equals total_scv (up to rounding)
    within_scv_alt  = within_scv_alt,    # should match within_scv
    total_var       = total_var,
    within_var      = within_var,
    between_var     = between_var,
    check_var       = within_var + between_var   # equals total_var
  )
  
  if (isTRUE(per_group)) {
    out$by_group <- data.frame(
      group = names(mu_g),
      n = as.integer(n_g),
      p = as.numeric(p_g),
      mu_g = as.numeric(mu_g),
      var_g = as.numeric(var_g),
      scv_g = as.numeric(scv_g),
      contrib_within_scv  = as.numeric(p_g * (mu_g / mu)^2 * scv_g),
      contrib_between_scv = as.numeric(p_g * ( (mu_g / mu) - 1 )^2)
    )
  }
  
  out
}


#######################################################################################################
# 2. , 3 [Credit: Domenico Moramarco]
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
############################################################################################################
# 4. 
# Pure vector-in / scalar-out BM between component
bm_between_vec <- function(y, g) {
  y <- as.numeric(y); g <- as.factor(g)
  N  <- length(y); if (N < 2) return(0)
  mu <- mean(y)
  mu_g <- tapply(y, g, mean)
  n_g  <- tapply(y, g, length)
  p    <- n_g / N
  sum(outer(p, p) * abs(outer(mu_g, mu_g, `-`))) / (2 * mu)
}

################################################################################################################
# 5. 
# --- main wrapper ---
# data: your data frame
# id_col, income_col, group_col: column names (strings)
# join_df: optional external df with per-id "BM" (between)
# join_id_col, join_between_col: column names in join_df
fleurbaey_pipeline <- function(
  data,
  id_col            = "id",
  income_col        = "Rating",
  group_col         = "Approval",
  join_df           = NULL,
  join_id_col       = "id",
  join_between_col  = "BM"       # <-- default now "BM"
) {
  id_sym   <- sym(id_col)
  inc_sym  <- sym(income_col)
  grp_sym  <- sym(group_col)

  safe_gwP  <- function(y, g) tryCatch(gini_within_add_path(y, g)[2], error = function(e) NA_real_)
  safe_gini <- function(y)    tryCatch(gini_coefficient(y),           error = function(e) NA_real_)

  # per-id: gini_wit_path and Res_new = G - GwP - GbA
  res_ids <- data %>%
    filter(is.finite(!!inc_sym), !is.na(!!grp_sym)) %>%
    group_by(!!id_sym) %>%
    group_modify(~{
      y <- .x[[income_col]]
      g <- .x[[group_col]]
      GwP <- safe_gwP(y, g)              # within (path)
      GbA <- bm_between_vec(y, g)        # between (add)
      G   <- safe_gini(y)                # overall Gini
      tibble(
        gini_wit_path = GwP,
        Res_new       = G - GwP - GbA    # R* = G - GwP - GbA
      )
    }) %>%
    ungroup() %>%
    rename(id = !!id_sym)

  if (is.null(join_df)) {
    return(list(res_ids = res_ids))
  }

  # Join with external per-id BM and compute WDP.Fleurbaey (using 0.5*Res_new as in your code)
  join_id_sym  <- sym(join_id_col)
  join_bm_sym  <- sym(join_between_col)

  fleurbaey <- join_df %>%
    select(!!join_id_sym, !!join_bm_sym) %>%
    rename(id = !!join_id_sym, Gini_between = !!join_bm_sym) %>%
    mutate(Gini_between = as.numeric(Gini_between)) %>%
    left_join(res_ids, by = "id") %>%
    mutate(WDP.Fleurbaey = ifelse(gini_wit_path < (Gini_between - 0.5 * Res_new), 1L, 0L))

  wdp_table   <- table(fleurbaey$WDP.Fleurbaey, useNA = "ifany")
  resnew_summ <- summary(fleurbaey$Res_new)

  list(
    res_ids      = res_ids,
    fleurbaey    = fleurbaey,
    wdp_table    = wdp_table,
    resnew_summ  = resnew_summ
  )
}

## apply the functions
# Just the per-id results (no join):
#out1 <- fleurbaey_pipeline(france_theil.df)
#out1$res_ids

# With your external table ic2res.df (must have columns id and Gini.between):
#out2 <- fleurbaey_pipeline(
#  data = france_theil.df,
#  join_df = ic2res.df,
#  id_col = "id",
#  income_col = "Rating",
#  group_col = "Approval",
#  join_id_col = "id",
#  join_between_col = "Gini.between"
#)

#out2$res_ids
#out2$fleurbaey
#out2$wdp_table
#out2$resnew_summ
####
## not run
##########################################################################################################
# 6. Calculate Decomposition
ic2decomp.fun <- function(i, data) {
  tryCatch({
    df <- data %>% 
      filter(id %in% i) %>% 
      mutate(Approval = factor(Approval)) 
    
    dec <- decompGEI(x = df$Rating, z = df$Approval, alpha = 1, ELMO = FALSE)
    th.total <- as.numeric(dec$ineq$index)
    th.within <- as.numeric(dec$decomp$within)
    th.between <- as.numeric(dec$decomp$between)
    # Compute Gini decomposition
    gin <- decompSGini(x = df$Rating, z = df$Approval, decomp = "YL", ELMO = FALSE)  # YL
    gi.total <- as.numeric(gin$ineq$index)
    gi.within <- as.numeric(gin$decomp$within)
    gi.between <- as.numeric(gin$decomp$between)
    BM <- decompSGini(x = df$Rating, z = df$Approval, decomp = "BM", ELMO = TRUE)  # YL
    # Compute Atkinson decomposition
    atk <- decompAtkinson(x = df$Rating, z = df$Approval, decomp = "BDA", epsilon = 1)
    atk.total <- as.numeric(atk$ineq$index)
    atk.within <- as.numeric(atk$decomp$within)
    atk.between <- as.numeric(atk$decomp$between)
    ## Squared Coefficient of Variation
    scv <- cv_decomp_bm(df, income = "Rating", group = "Approval", per_group = TRUE)
    scv.total <- scv$total_scv
    scv.within <- scv$within_scv
    scv.between <- scv$between_scv 
    ## Decomposition of the variance
 #   variance <- var_decomp(df, income="Rating", group="Approval",
 #                 type="sample",per_group=TRUE)
 #   var.total <- variance$total_var
 #   var.within <- variance$within_var
 #   var.between <- variance$between_var
    ## Mean-Log Deviation
 #   mld <- dineq::mld_decomp(x = df$Rating, z = df$Approval)
 #   mld.total <- mld$mld_decomp$mld_total
 #   mld.within <- mld$mld_decomp$mld_within
 #   mld.between <- mld$mld_decomp$mld_between
    # Store results in a data frame
    res.df <- data.frame(
      id = unique(df$id),
      Theil.tot = th.total,
      Theil.within = th.within,
      Theil.between = th.between,
      Gini.tot = gi.total,
      Gini.within = gi.within,
      Gini.between = gi.between,
      Atkinson.tot = atk.total,
      BM.between = as.numeric(BM$decomp$between), ## Gini decompositino with BM
      Atkinson.within = atk.within,
      Atkinson.between = atk.between,
  #    MLD.tot = mld.total,
  #    MLD.within = mld.within,
  #    MLD.between = mld.between
      SCV.tot = scv.total,
      SCV.within = scv.within,
      SCV.between = scv.between#,
#      VAR.total = var.total,
#      VAR.within = var.within,
#      VAR.between = var.between
    )
    return(res.df)
  }, error = function(e) {
    message(paste("Skipping id:", i, "due to error:", e$message))
  })
}

save.image("~/Documents/Research/Dichotomous/github/Dichotomous/AuxFunctions.RData")
