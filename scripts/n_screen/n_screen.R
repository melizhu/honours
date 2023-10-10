args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) { # no arguments provided
  stop("Parameter for batch number needs to be provided.")
} else {
  batch <- as.integer(args[1])
}

source("scripts/NEW_ALL_functions.R")

# Define the number of different kinds of mutations
n <- 2^(0:7)

# parameters:
K <- 10^6
A <- 1
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

for (i in 1:length(n)) {
  cat(paste0("Running simulations for n=", n[i], ".\n"))
  zeta_m <- rgamma(n[i], shape = zeta_shape, scale = zeta_scale) + 1
  #The code with building data frame for S and all the mutations
  strains <- data.frame(name = c("S", paste0("M", 1:n[i])),
                        psi = psi,
                        phi = phi,
                        kappa = kappa,
                        zeta = c(1, zeta_m))
  sims[[i]] <- list()
  sims[[i]]$strains <- strains
  sims[[i]]$generalinput <- c(K = K,
                              A = A,
                              miu = miu,
                              Sini = Sini,
                              zeta_shape = zeta_shape,
                              zeta_scale = zeta_scale,
                              rseed = rseed)

  sims[[i]]$results <- replicate(m, simulate(strains, K = K, A = A, miu = miu, Sini = K))
}

save(sims, file = paste0("output/n_screen/SimResults_n_screen_rep_", batch, ".RData"))
