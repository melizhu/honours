#orginal function

#change the growth rate of S to test plot 4 of sims 1 in t1 of sim=3
#change the death rate
library(adaptivetau)
library(vegan)

# simulate is used for generating simulations 
simulate <- function(strains, K, A, miu, tf=1000, Sini)
{
  n <- nrow(strains)
  states<-strains$name
  #m is the death rate function
  #m is the reduction in bacteria fitness according to antibiotic concentration A
  #zeta is MIC(minimal inhibitory concentration) of the drug
  #psi is the growth rate in the absences of antibiotics 
  #phi is the maximum reduction in the fitness
  #kappa is the slope
  m <- function(A, phi, psi, zeta, kappa) {
    phi * ((A / zeta) ^ kappa) / ((A / zeta) ^ kappa - (psi - phi) / psi)
  }
  
  #built the transition rate at loop1 as events
  
  transitions <- list()
  for (i in 1){ # build S growth and dies
    transitions[[i]] <- setNames(rep(0, length(states)), states)
    transitions[[i]][1]<-1
    transitions[[2*i]] <- setNames(rep(0, length(states)), states)
    transitions[[2*i]][1]<--1
  }
  for (i in 3:(n+1)) { # Define S to mutated to M1 M2 M3 etc...
    transitions[[i]] <- setNames(rep(0, length(states)), states)
    transitions[[i]][1]<--1
    transitions[[i]][i-1]<-1
    
  }
  
  for (i in 2:n){
    ##Define M1,M2,M3 etc growth
    transitions[[i+n]] <- setNames(rep(0, length(states)), states)
    transitions[[i+n]][i]<-1
    #Define M1,M2,M3 etc dies
    transitions[[i+2*n-1]] <- setNames(rep(0, length(states)), states)
    transitions[[i+2*n-1]][i]<--1
  }
  
  # Function to calculate transition rates, given variables and parameters
  rates <- function(x, pars, t) {
    
    with(as.list(c(pars$general, pars$S)), {
      N <- sum(unname(x[strains$name]))
      #N define the whole population
      Sdie <<- x["S"]*m(A, phi, psi, zeta, kappa)# s dies
      Sgrowth <<- x["S"]*(psi*(1-(N)/K)*(1-miu)) # s grows
      if(Sgrowth < 0) { # to avoid the condition that we had negative growth 
        Sdie <<- Sdie - Sgrowth # add to death population
        Sgrowth <<- 0  # set the negative growth to zero
       
        
      }
      #set sgrowth <Sdie? if sgrowth=Sdie
      #if (Sgrowth >Sdie){
        #Sgrowth<<-0
      #}
    })
    
    StoM<- numeric(n-1) # susceptible population mutated to mutation population 
    growthM <- numeric(n-1) # mutation grows 
    dieM <- numeric(n-1) # mutation dies
    
    for (i in 1:(n-1)) {
      with(as.list(c(states, pars$general,pars$S)), {
        N <- sum(unname(x[strains$name]))
        StoM[i] <<- x["S"]*(psi*(1-(N)/K)*(miu/(n-1))) # susceptible population mutated to mutation population 
        if(StoM[i] < 0) { # to avoid the condition that we had negative population mutated 
          StoM[i]<<-0 # set the negative to zero
        }
      })
    }
    for (i in 1:(n-1)) {
      with(as.list(c(states, pars$general, pars[[paste0("M", i)]])), {
        dieM[i] <<- x[paste0("M", i)]* m(A, phi, psi, zeta, kappa)
        N <- sum(unname(x[strains$name]))
        growthM[i] <<- x[paste0("M", i)]*psi*(1-N/K) # mutation grows
        if(growthM[i] < 0) { # to avoid the condition that we had negative growth 
          dieM[i] <<- dieM[i] - growthM[i] # add to death population
          growthM[i] <<- 0 # set the negative growth to zero
        }
      })
    }
    
    if(any(growthM < 0)) { # check to avoid the condition that we had negative growth 
      print(t)
      print(x)
      print(growthM)
    }
    
    return(c(Sgrowth,
             Sdie,
             StoM,
             growthM,
             dieM ))
    
  }
  
  
  pars <- list(general = list(K = K, A = A, miu =miu))
  
  for (i in 1:n) {
    pars[strains$name[i]] <- list(list(psi = strains$psi[i],
                                       phi = strains$phi[i],
                                       kappa =strains$kappa[i] ,
                                       zeta = strains$zeta[i]))
  }
  
  
  initial <- c(S=Sini) # initialize S
  
  
  for (i in 1:(n-1)) {
    initial[paste0("M", i)] <- 0  # initialize mutations
  }
  
  r <- ssa.adaptivetau(initial, transitions, rates, pars, tf)
  
  return(r)
}

#test the function with input of n=5 sims[[3]] results[[1:6]]
# Define the number of different kinds of mutations
n <- 4

# parameters:
K <- 10^6 # carrying capacity
A <- 1 # antibiotic concentration #0
miu <- 10^(-7) # mutation rate 
Sini <- K # initial population

psi <- 0.7 # growth rate in the absence of antibiotics 
phi <- 5 # maximum reduction in fitness
kappa <- 2.5 # slope
zeta_shape <- 2 # random gamma distribution shape for generating mutation MIC
zeta_scale <- 5 # random gamma distribution scale for generating mutation MIC

# Initialize a list to store the data frames
sim_test <- list()
# number of simulations
m <- 10

# make a random seed:
rseed <- 9.3873e+05
set.seed(rseed)

# for generating the simulations

  cat(paste0("Running simulations for n=", n, ".\n"))
  #the zeta values for mutations are random gamma distribution and greater than 1
  zeta_m <- rgamma(n, shape = zeta_shape, scale = zeta_scale) + 1
  #The code with building data frame for S and all the mutations
  strains <- data.frame(name = c("S", paste0("M", 1:n)),
                        psi = psi,
                        phi = phi,
                        kappa = kappa,
                        zeta = c(1, zeta_m))
  #build a list for saving simulation inputs and results called sims
  sim_test <- list()
  #sims strains contain all the PD parameters such as psi, phi, kappa and zeta
  sim_test$strains <- strains
  #sims generalinput contain all the other variables 
  sim_test$generalinput <- c(K = K, 
                              A = A,
                              miu = miu,
                              Sini = Sini,
                              zeta_shape = zeta_shape,
                              zeta_scale = zeta_scale,
                              rseed = rseed)
  #sims results contain the simulation results 
  sim_test$results <- replicate(m, simulate(strains, K = K, A = A, miu = miu, Sini = K))

#plotting the results 
library(tidyverse)
#here we pick sim[[3]] results[[4]]
#i <- 3
j <- 3
#plot this particular simulation result
simtestplot <- sim_test$results[[j]] |>
  as_tibble() |>
  pivot_longer(-time)
ggplot(simtestplot) +
  geom_line(aes(x = time, y = log(value), colour = name))