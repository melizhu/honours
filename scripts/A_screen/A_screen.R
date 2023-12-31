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
K <- 10^6
A <-seq(0.8,1.3,by=0.05)# 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 [1:11]
miu <- 10^(-7)
Sini <- K

psi <- 0.7
phi <- 5.0
kappa <- 2.5
zeta_shape <- 2
zeta_scale <- 5

# Initialize a list to store the data frames
sims <- list()
# number of simulations
m <- 1000

# make a random seed:
rseed <- sample(1:1e6, 1)
set.seed(rseed)

for (i in 1:length(A)) {
  cat(paste0("Running simulations for A=", A[i], ".\n"))
  zeta_m <- rgamma(n, shape = zeta_shape, scale = zeta_scale) + 1
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
                              A = A[i],
                              miu = miu,
                              Sini = Sini,
                              zeta_shape = zeta_shape,
                              zeta_scale = zeta_scale,
                              rseed = rseed)
  #sims results contain the simulation results
  sims[[i]]$results <- replicate(m, simulate(strains, K = K, A = A[i], miu = miu, Sini = Sini))
}
#save the simulation results to a file
save(sims, file = paste0("output/A_screen/SimResults_A_screen_rep_", batch, ".RData"))
