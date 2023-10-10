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
K <- as.integer(10^seq(5,9,by=0.25)) # K->100000     177827     316227     562341    1000000    1778279    3162277    5623413   10000000
# 17782794   31622776   56234132  100000000  177827941  316227766  562341325 1000000000 [1:17]

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

for (i in 1:length(K)) {
  cat(paste0("Running simulations for K=", K[i], ".\n"))
  zeta_m <- rgamma(n, shape = zeta_shape, scale = zeta_scale) + 1
  #The code with building data frame for S and all the mutations
  strains <- data.frame(name = c("S", paste0("M", 1:n)),
                        psi = psi,
                        phi = phi,
                        kappa = kappa,
                        zeta = c(1, zeta_m))
  sims[[i]] <- list()
  sims[[i]]$strains <- strains
  sims[[i]]$generalinput <- c(K = K[i],
                              A = A,
                              miu = miu,
                              Sini = Sini[i],
                              zeta_shape = zeta_shape,
                              zeta_scale = zeta_scale,
                              rseed = rseed)
  
  sims[[i]]$results <- replicate(m, simulate(strains, K = K[i], A = A, miu = miu, Sini = Sini[i]))
}

save(sims, file = paste0("output/p_screen/SimResults_pop_screen_rep_", batch, ".RData"))
