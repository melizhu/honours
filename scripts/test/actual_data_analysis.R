#Here is a example code that I used to do analyzing the data
#Here I was doing analyzing for p_screen to find what are those median values
#S here has the same value initial population size I choose
S <-as.integer(10^seq(5, 9, by = 0.25))
S<- as.numeric(S)
# Create an empty vector to store the median values
save_median <- numeric(length(S))

# Loop through the values in 'n'
for (i in 1:length(S)) {
  # Find the first row number where summary$n == n[i]
  rownumber1 <- min(which(summary$K == S[i]))
  
  # Find the last row number where summary$n == n[i]
  rownumber2 <- max(which(summary$K == S[i]))
  
  # Subset the data frame to rows where summary$n == n[i]
  subset_data <- summary[rownumber1:rownumber2, ]
  
  # Calculate the median of the 'ed_average_without_wildtype' column for this subset
  save_median[i] <- median(subset_data$sims_with_mutations)
}

# 'save_median' now contains the median values for each subset 
save_median