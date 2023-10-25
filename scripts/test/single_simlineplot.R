#plot single line plots of simulation
library(tidyverse)
#here we pick sim[[3]] results[[1]]
i <- 3
j <- 1
#plot this particular simulation result
singlesim <- sims[[i]]$results[[j]] |>
  as.tibble() |>
  pivot_longer(-time)
ggplot(singlesim) +
  geom_line(aes(x = time, y = value, colour = name))

#create a new big table for all the simulations
t1<-bind_rows(lapply(sims, function(sim){
  bind_rows(lapply(sim$results, function(results){
    as.tibble(results) |> pivot_longer(-time, names_to = "strain", values_to = "pop")
  }), .id="replicate")
}), .id="sim")

# plot line plots for multiple simulations
t1 %>%
  # here we plot 1 to 6 simulations from total 1000 simulations
  filter(sim==3, replicate %in% 1:6) %>%
  ggplot() +
  geom_line(aes(x = time, y = log(pop), colour = strain)) +
  labs(x = "time", y = "ln(population size)") +
  facet_wrap(~replicate, ncol = 3)+
  theme_bw()


# Extract the number of columns from the first results data frame
num_cols <- ncol(sims[[i]]$results[[1]])
# Create the freq_T matrix with the appropriate number of columns
freq_T <- data.frame(matrix(0, nrow = 1000, ncol = num_cols))
colnames(freq_T) <- colnames(sims[[i]]$results[[1]])

for (j in 1:length(sims[[i]]$results)) {
  #entropy and ed_average with S
  final_numbers <- tail(sims[[i]]$results[[j]][, -1], 1)
  
  
  #entropy fill in freq_T with the winner=1
  strain_name <- colnames(final_numbers)[which.max(final_numbers)]
  freq_T[j , which(colnames(freq_T) == strain_name)] <- 1
  
}

#generate the "winning" distribution for each "winner" 
distribution<- apply(freq_T[, 2:ncol(freq_T)], 2, function(x) sum(x)/ nrow(freq_T))

# Plotting a histogram
barplot(distribution, col = "skyblue", main = "Winning distribution of 1000 simulations",
        xlab = "Categories", ylab = "Values")


