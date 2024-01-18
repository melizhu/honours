args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) { # no arguments provided
  stop("Parameter for batch number needs to be provided.")
} else {
  batch <- as.integer(args[1])
}

source("scripts/NEW_ALL_functions.R")

# Define the number of different kinds of mutations
n <- 32

# parameters:
K <- 10^6 # carrying capacity
A <- 1 # antibiotic concentration
miu <- 10^(-7) # mutation rate 
Sini <- K # initial population

psi <- 0.7 # growth rate in the absence of antibiotics 
phi <- 5.0  # maximum reduction in fitness
kappa <- 2.5 # slope

#zeta screen
zeta_shape <- 2^seq(-6, 9)#c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)# random gamma distribution shape for generating mutation MIC
zeta_scale <- 10/(zeta_shape) # random gamma distribution scale for generating mutation MIC

# Initialize a list to store the data frames
sims <- list()
zeta_save<-list()
# number of simulations
m <- 1000

# make a random seed:
rseed <- sample(1:1e6, 1)
set.seed(rseed)

#write loop for generating the simulations
for (i in 1:length(zeta_shape)) {
  cat(paste0("Running simulations for MIC(zeta_shape)=", zeta_shape[i], ".\n"))
  #the zeta values for mutations are random gamma distribution and greater than 1
  zeta_m<- rgamma(n, shape = zeta_shape[i], scale = zeta_scale[i]) + 1
  #save zeta_m each time
  zeta_save[[i]]<-zeta_m
  #The code with building data frame for S and all the mutations
  strains <- data.frame(name = c("S", paste0("M", 1:n)),
                        psi = psi,
                        phi = phi,
                        kappa = kappa,
                        zeta = c(1, zeta_m))
  #build a list for saving simulation inputs and results called sims
  sims[[i]] <- list()
  #sims strains contain all the PD parameters such as psi, phi, kappa and zeta
  sims[[i]]$strains <- strains
  #sims generalinput contain all the other variables 
  sims[[i]]$generalinput <- c(K = K, 
                              A = A,
                              miu = miu,
                              Sini = Sini,
                              zeta_shape = zeta_shape[i],
                              zeta_scale = zeta_scale[i],
                              rseed = rseed)
  #sims results contain the simulation results 
  sims[[i]]$results <- replicate(m, simulate(strains, K = K, A = A, miu = miu, Sini = K))
}

#save the simulation results to a file 
save(sims, file = paste0("output/z_screen/SimResults_z_screen_rep_", batch, ".RData"))

# checking the PD curves:
# plot_PDs(sims[[1]]$strains)
