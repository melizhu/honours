
source("scripts/NEW_ALL_functions.R")

#set n is the same for the HR and LR
# Define the number of different kinds of mutations
n <- 4


#choose small n big s0 and high A low zeta shape for high repeatability
# parameters:
K <- c(10^9,10^6) # carrying capacity
A <-c(1.3,1) # antibiotic concentration
miu <- 10^(-7) # mutation rate 
Sini <- K # initial population

psi <- 0.7 # growth rate in the absence of antibiotics 
phi <- 5.0  # maximum reduction in fitness
kappa <- 2.5 # slope
zeta_shape <- c(2^(-6),2) # random gamma distribution shape for generating mutation MIC
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

#plot the sim_1 results [2][3][4] with x-axis is time and y axis is the population 

#stick those two plot together and put the HR on the left and LR on the right with title 
for (i in 1:2) {
  for(j in 2:4) {
    singlesim <- sims_1[[i]]$results[[j]] %>%
      as_tibble() %>%
      pivot_longer(-time)
    
    assign(paste0("plot", (i - 1) * 3 + j - 1), # To access the plots, use plot1, plot2, ..., plot6
           ggplot(singlesim) +
             geom_line(aes(x = time, y = log(value), colour = name)) +
             labs(title = paste("Simulation", (i - 1) * 3 + j - 1)))
  }
}


# Arrange plots in a grid
example_figure<-grid.arrange(
  arrangeGrob(
    plot1, plot2, plot3,
    nrow = 3, ncol = 1    #and they are in 1:3 1 for colomn and 3 for row in the figure
  ),
  arrangeGrob(
    plot4, plot5, plot6,
    nrow = 3, ncol = 1  #and they are in 1:3 1 for colomn and 3 for row in the figure
  ),
  top = textGrob("HR", x = 0.25, y = 0.34, just = "center", gp = gpar(fontsize = 20)),
  nrow = 1
)

# Add LR title manually
plots<-arrangeGrob(
  example_figure,
  grid.text("LR", x = 0.75, y = 0.98, just = "center", gp = gpar(fontsize = 20))
)
  
  
ggsave("plots/example_figure.pdf", plots, width = 3, height = 3)