library(tidyverse)
library(adaptivetau)
library(vegan)
#reset up all the useful parameters just in case it is other values
l <-1
m <- 1000
repfile<-32

zeta_shape <- 2^seq(-6, 8)
zeta_scale <- 10/(zeta_shape)

# Create empty data frame for both entropy values , average_ed, Bray-Curtis_distance and rate of mutations in sims values
summary <- data.frame(zeta_shape = rep(zeta_shape, each = repfile),
                      PDtable = rep(1:repfile, length(zeta_shape)),
                      entropy = NA,
                      entropy_without_wildtype=NA,
                      ed_average=NA,
                      ed_average_without_wildtype=NA,
                      bray_distance_average=NA, #Bray-Curtis_distance average
                      bray_average_without_wildtype=NA, #Bray-Curtis_distance average _without_wildtype
                      sims_with_mutations=NA)






analyseSims_zeta<-function(sims){
  entropy_function <- function(p) { # entropy function
    p <- p[p > 1e-10]
    - sum(p * log(p, base = 2))}
  
  
  # Extract the number of columns from the first results data frame
  num_cols <- ncol(sims$results[[1]])
  # Create the freq_T matrix with the appropriate number of columns
  freq_T <- data.frame(matrix(0, nrow = m, ncol = num_cols))
  colnames(freq_T) <- colnames(sims$results[[1]])
  # Create the eucli and bray matrix with the appropriate number of columns 
  eucli_old<- data.frame(matrix(0, nrow = m, ncol = num_cols-1))
  # Create the eucli and bray matrix without S with the appropriate number of columns 
  eucli_new<- data.frame(matrix(0, nrow = m, ncol = num_cols-2))
  
  
  # build a empty house for the mutation_grwoth_count 
  mutation_grwoth_count<-vector("numeric",length=m)
  # build mutations_grwoth
  mutations_grwoth<-0
  
  # Fill freq_T matrix
  for (j in 1:length(sims$results)) {
    #entropy and ed_average with S
    final_numbers <- tail(sims$results[[j]][, -1], 1)
    
    
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
  sims_summary$entropy<- entropy_function(p)
  
  # calculate entropy_without_wildtype
  q <- c(unname(freq_new))
  sims_summary$entropy_without_wildtype <- entropy_function(q)
  
  
  
  #calculate the ed_distance_average
  # Remove rows with NA values
  eucli_old_cleaned <- na.omit(eucli_old)
  ed_distance<-vegdist(eucli_old_cleaned,method="euclidean") #euclidean
  average_ed <- mean(as.vector(ed_distance), na.rm = TRUE) 
  sims_summary$ed_average <-  average_ed
  
  #calculate the ed_without_wildtype_average
  # Remove rows with NA values 
  eucli_new_cleaned <- na.omit(eucli_new)
  distance_without_wildtype<-vegdist(eucli_new_cleaned,method="euclidean") 
  average_without_wildtype <- mean(as.vector(distance_without_wildtype), na.rm = TRUE)
  sims_summary$ed_average_without_wildtype <-  average_without_wildtype
  
  
  #calculate the Bray-Curtis_distance average
  # Remove rows with NA values
  bray_distance<-vegdist(eucli_old_cleaned,method="bray") 
  average_bray_distance <- mean(as.vector(bray_distance), na.rm = TRUE) 
  sims_summary$bray_distance_average <-  average_bray_distance
  
  #calculate the bray_distancewithout_wildtype
  # Remove rows with NA values
  bray_distance_without_wildtype<-vegdist(eucli_new_cleaned,method="bray") #euclidean
  average_bray_without_wildtype <- mean(as.vector( bray_distance_without_wildtype), na.rm = TRUE) 
  sims_summary$bray_average_without_wildtype <-  average_bray_without_wildtype
  
  #calculate the sims_with_mutations
  mutations_grwoth<-sum(mutation_grwoth_count)/m
  sims_summary$sims_with_mutations <-mutations_grwoth
  
  
  return(sims_summary)
}
start.time <- Sys.time()

#go through all the sims from rep1 
for (i in 1:32) {
  for(j in 1:15){
  # Load the file
 load(paste0("output/z_screen/SimResults_z_screen_rep_",i,"_zeta_", j, ".RData"))
  # Create empty data frame
  sims_summary <- data.frame(zeta_shape = zeta_shape[i],
                             PDtable = NA,
                             entropy = NA,
                             entropy_without_wildtype=NA,
                             ed_average=NA,
                             ed_average_without_wildtype=NA,
                             bray_distance_average=NA,
                             bray_average_without_wildtype=NA,
                             sims_with_mutations=NA)
  # Analyze the loaded file using the analyseSims function
  analysis_result <- analyseSims_zeta(sims)  
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i & summary$zeta_shape == zeta_shape[j], "entropy"] <- analysis_result$entropy
  # Fill analysis_result$ into the summary that has PDtable value equal to i
  summary[summary$PDtable == i & summary$zeta_shape == zeta_shape[j], "entropy_without_wildtype"] <- analysis_result$entropy_without_wildtype
  # Fill analysis_result$ into the summary that has PDtable value equal to i
  summary[summary$PDtable == i & summary$zeta_shape == zeta_shape[j], "ed_average"] <- analysis_result$ed_average
  # Fill analysis_result$ into the summary that has PDtable value equal to i
  summary[summary$PDtable == i & summary$zeta_shape == zeta_shape[j], "ed_average_without_wildtype"] <- analysis_result$ed_average_without_wildtype
  # Fill analysis_result$ into the summary that has PDtable value equal to i
  summary[summary$PDtable == i & summary$zeta_shape == zeta_shape[j], "bray_distance_average"] <- analysis_result$bray_distance_average
  # Fill analysis_result$ into the summary that has PDtable value equal to i
  summary[summary$PDtable == i & summary$zeta_shape == zeta_shape[j], "bray_average_without_wildtype"] <- analysis_result$bray_average_without_wildtype
  # Fill analysis_result$ into the summary that has PDtable value equal to i
  summary[summary$PDtable == i & summary$zeta_shape == zeta_shape[j], "sims_with_mutations"] <- analysis_result$sims_with_mutations
  }
}


end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken


#plot entropy
plot1<-ggplot(summary, aes(x = as.factor(zeta_shape), y = entropy)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "z(zeta_shape)", y = "Entropy") +
  theme_bw()

#plot ed_average
plot2<-ggplot(summary, aes(x = as.factor(zeta_shape), y = ed_average)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "z(zeta_shape)", y = "Eucliean distance") +
  theme_bw()
#plot entropy_without_wildtype
plot3<-ggplot(summary, aes(x = as.factor(zeta_shape), y = entropy_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "z(zeta_shape)", y = "Entropy") +
  theme_bw()
#plot ed_average_without_wildtype
plot4<-ggplot(summary, aes(x = as.factor(zeta_shape), y = ed_average_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "z(zeta_shape)", y = "Euclidean Distance") +
  theme_bw()

#plot bray_average
plot5<-ggplot(summary, aes(x = as.factor(zeta_shape), y = bray_distance_average)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "z(zeta_shape)", y = "Bray-Curtis Distance") +
  theme_bw()
#plot bray_average_without_wildtype
plot6<-ggplot(summary, aes(x = as.factor(zeta_shape), y = bray_average_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "z(zeta_shape)", y = "Bray-Curtis Distance ") +
  theme_bw()
#plot mutation sims_with_mutations  frequency 
plot7<-ggplot(summary, aes(x = as.factor(zeta_shape), y = sims_with_mutations )) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "z(zeta_shape)", y = "P(mutation spread)") +
  theme_bw()
#save plots
ggsave("plots/z_screen/entropy_zscreen.pdf", plot1, width = 8, height = 8)
ggsave("plots/z_screen/ed_average_zscreen.pdf", plot2, width = 8, height = 8)
ggsave("plots/z_screen/entropy_without_wildtype_zscreen.pdf", plot3, width = 8, height = 8)
ggsave("plots/z_screen/ed_average_without_wildtype_zscreen.pdf", plot4, width = 8, height = 8)
ggsave("plots/z_screen/bray_average_zscreen.pdf", plot5, width = 8, height = 8)
ggsave("plots/z_screen/bray_average_without_wildtype_zscreen.pdf", plot6, width = 8, height = 8)
ggsave("plots/z_screen/sims_with_mutations_zscreen.pdf", plot7, width = 8, height = 8)

# save the database
save(summary, file = "output/z_screen/summary_z.RData")

#combine plots 
library(cowplot)
#Plot all graphs contains mutants
plot_all_mutants<-plot_grid(plot7,plot3,plot4,plot6, labels = c("(a)", "(b)","(c)", "(d)"), hjust=0, vjust=1.3)
#Plot all graphs contains S and mutants
plot_all_S<-plot_grid(plot1,plot2,plot5, labels = c("(a)", "(b)","(c)"), hjust=0, vjust=1.3)

#save plots of the mutants
ggsave("plots/z_screen/plot_all_mutants.pdf", plot_all_mutants, width = 18, height = 15)
#save plots of S and mutants
ggsave("plots/z_screen/plot_all_S.pdf", plot_all_S, width = 18, height=15)


