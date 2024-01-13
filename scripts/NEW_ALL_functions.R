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

  r <- ssa.adaptivetau(initial, transitions, rates, pars, tf, tl.params = list(epsilon=0.01))

  return(r)
}

# analyseSims is used for analyzing the simulation results and generating results 
analyseSims<-function(sims){
  entropy_function <- function(p) { # entropy function
    p <- p[p > 1e-10]
    - sum(p * log(p, base = 2))}

  
  for (i in 1:length(sims)) {
    
    # Extract the number of columns from the first results data frame
    num_cols <- ncol(sims[[i]]$results[[1]])
    # Create the freq_T matrix with the appropriate number of columns
    freq_T <- data.frame(matrix(0, nrow = m, ncol = num_cols))
    colnames(freq_T) <- colnames(sims[[i]]$results[[1]])
    # Create the eucli and bray matrix with the appropriate number of columns 
    eucli_old<- data.frame(matrix(0, nrow = m, ncol = num_cols-1))
    # Create the eucli and bray matrix without S with the appropriate number of columns 
    eucli_new<- data.frame(matrix(0, nrow = m, ncol = num_cols-2))
    
    
    # build a empty house for the mutation_grwoth_count 
    mutation_grwoth_count<-vector("numeric",length=m)
    # build mutations_grwoth
    mutations_grwoth<-0
    
    # Fill freq_T matrix
    for (j in 1:length(sims[[i]]$results)) {
      #entropy and ed_average with S
      final_numbers <- tail(sims[[i]]$results[[j]][, -1], 1)
      
    
      #entropy fill in freq_T with the winner=1
      strain_name <- colnames(final_numbers)[which.max(final_numbers)]
      freq_T[j , which(colnames(freq_T) == strain_name)] <- 1
      
      
      
      #ed_average with S fill in the matrix eucli_old
      total_old <- sum(final_numbers)
      if (total_old == 0) {
        freq_old <- NA
      } else {
        freq_old <-  final_numbers / total_old
      }
      
      
      eucli_old[j,]<-freq_old
      colnames(eucli_old) <- colnames(final_numbers)
      
      #ed_average without S fill in the matrix eucli_new
      final_number_without_wildtype<-final_numbers[, -1]
      
      total_new <- sum(final_number_without_wildtype)
      if (total_new == 0) {
        freq_new <- NA
      } else {
        freq_new <-  final_number_without_wildtype / total_new
      }
      
      eucli_new[j,]<-freq_new
      colnames(eucli_new) <- colnames(final_number_without_wildtype)
      
      #fill in the mutation_grwoth_count
      if (total_new >0) {
        mutation_grwoth_count[j]<-1
      } else{
        mutation_grwoth_count[j]<-0
      }
    }
    
    # Find rows with S = 1 and exclude it to create freq_T_new
    freq_T_new <- freq_T[freq_T$S != 1, ]
    # Calculate frequencies
    freq<- apply(freq_T[, 2:ncol(freq_T)], 2, function(x) sum(x) / nrow(freq_T))
    # Calculate frequencies for freq_T_new
    freq_new <- apply(freq_T_new[, 2:ncol(freq_T_new)], 2, function(x) sum(x) / nrow(freq_T_new))
    
    # Calculate entropy
    p <- c(unname(freq))
    sims_summary$entropy[i] <- entropy_function(p)
    
    # calculate entropy_without_wildtype
    q <- c(unname(freq_new))
    sims_summary$entropy_without_wildtype[i] <- entropy_function(q)
    
    
    
    #calculate the ed_distance_average
    # Remove rows with NA values
    eucli_old_cleaned <- na.omit(eucli_old)
    ed_distance<-vegdist(eucli_old_cleaned,method="euclidean") #euclidean
    average_ed <- mean(as.vector(ed_distance), na.rm = TRUE) 
    sims_summary$ed_average[i] <-  average_ed
    
    #calculate the ed_without_wildtype_average
    # Remove rows with NA values 
    eucli_new_cleaned <- na.omit(eucli_new)
    distance_without_wildtype<-vegdist(eucli_new_cleaned,method="euclidean") 
    average_without_wildtype <- mean(as.vector(distance_without_wildtype), na.rm = TRUE)
    sims_summary$ed_average_without_wildtype[i] <-  average_without_wildtype
    
    
    #calculate the Bray-Curtis_distance average
    # Remove rows with NA values
    bray_distance<-vegdist(eucli_old_cleaned,method="bray") 
    average_bray_distance <- mean(as.vector(bray_distance), na.rm = TRUE) 
    sims_summary$bray_distance_average[i] <-  average_bray_distance
    
    #calculate the bray_distancewithout_wildtype
    # Remove rows with NA values
    bray_distance_without_wildtype<-vegdist(eucli_new_cleaned,method="bray") #euclidean
    average_bray_without_wildtype <- mean(as.vector( bray_distance_without_wildtype), na.rm = TRUE) 
    sims_summary$bray_average_without_wildtype[i] <-  average_bray_without_wildtype
    
    #calculate the sims_with_mutations
    mutations_grwoth<-sum(mutation_grwoth_count)/m
    sims_summary$sims_with_mutations[i] <-mutations_grwoth
    
  }
  return(sims_summary)
}