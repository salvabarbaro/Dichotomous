library(ineq)
library(dineq)
library(haven)
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
#     weighted_gini_add <- weight_g_add      + weight_g_add * gini_g
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
  
