library(ggplot2)
library(gridExtra)
library(tidyverse)
library(grid)
library(dplyr)

source("scripts/NEW_ALL_functions.R")

#set n is the same for the HR and LR
# Define the number of different kinds of mutations
n <- 4


#choose small n big s0 and high A low zeta shape for high repeatability
# parameters:
K <- c(10^9,10^6) # carrying capacity
A <-c(0.7,0.7) # antibiotic concentration
miu <- 10^(-7) # mutation rate 
Sini <- K # initial population

psi <- 0.7 # growth rate in the absence of antibiotics 
phi <- 5.0  # maximum reduction in fitness
kappa <- 2.5 # slope
zeta_shape <- c(2,2) # random gamma distribution shape for generating mutation MIC
zeta_scale <- 10/(zeta_shape) # random gamma distribution scale for generating mutation MIC

# Initialize a list to store the data frames
sims_1 <- list()
# number of simulations
m <- 10



#load all the sim_1 data by using the NEW_ALL_functions.R
# make a random seed:
rseed <- sample(1:1e6, 1)
set.seed(rseed)

for (i in 1:2) {
  #the zeta values for mutations are random gamma distribution and greater than 1
  zeta_m <- rgamma(n, shape = zeta_shape[i], scale = zeta_scale[i]) + 1
  #The code with building data frame for S and all the mutations
  strains <- data.frame(name = c("S", paste0("M", 1:n)),
                        psi = psi,
                        phi = phi,
                        kappa = kappa,
                        zeta = c(1, zeta_m))
  #build a list for saving simulation inputs and results called sims
  sims_1[[i]] <- list()
  #sims strains contain all the PD parameters such as psi, phi, kappa and zeta
  sims_1[[i]]$strains <- strains
  #sims generalinput contain all the other variables 
  sims_1[[i]]$generalinput <- c(K = K[i], 
                              A = A[i],
                              miu = miu,
                              Sini = Sini[i],
                              zeta_shape = zeta_shape[i],
                              zeta_scale = zeta_scale[i],
                              rseed = rseed)
  #sims results contain the simulation results 
  sims_1[[i]]$results <- replicate(m, simulate(strains, K = K[i], A = A[i], miu = miu, Sini = K[i]))


}


# Initialize an empty list to store the results
result_list <- list()

# Loop over the indices and replicate numbers
for (i in 1:2) {
  for (j in 1:3) {
    # Extract the data
    singlesim <- sims_1[[i]]$results[[j]] %>%
      as_tibble() %>%
      pivot_longer(-time) %>%
      filter(time <= 300)
    
    # Create the data frame
    replicate <- rep(paste("Replicate", j))
    parameter_set <- rep(ifelse(i == 1, "High Repeatability", "Low Repeatability"))
    singlesim_new <- data.frame(singlesim = singlesim, replicate = replicate, parameter_set = parameter_set)
    
    # Store the result in the list
    result_list[[length(result_list) + 1]] <- singlesim_new
  }
}

# Combine the results into a single data frame
S_together <- do.call(rbind, result_list)


ggplot(S_together) +
  geom_line(aes(x = singlesim.time, y = singlesim.value, colour = singlesim.name), show.legend = FALSE) +
  scale_y_log10()+
  scale_color_manual(values = c("purple", "blue", "red", "orange","black")) +
  theme_bw()+
  ylab("population")+
  xlab("time")+
  facet_grid(vars(replicate),vars(parameter_set))
  
ggsave("example_figure_new.jpg",example_figure, dpi = 300,width = 20, height =15 )
