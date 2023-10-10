library(tidyverse)

#reset up all the useful parameters just in case it is other values
l <-1
m <- 1000
repfile<-32

#K is the changing parameter in the pop_screen #change this make it looks nice 
A <-seq(0.8,1.3,by=0.05)


# Load the necessary functions
source("scripts/NEW_ALL_functions.R")
# Create empty data frame for both entropy values , average_ed, Bray-Curtis_distance and rate of mutations in sims values
summary <- data.frame(A = rep(A, each = repfile),
                      PDtable = rep(1:repfile, length(A)),
                      entropy = NA,
                      entropy_without_wildtype=NA,
                      ed_average=NA,
                      ed_average_without_wildtype=NA,
                      bray_distance_average=NA, #Bray-Curtis_distance average
                      bray_average_without_wildtype=NA, #Bray-Curtis_distance average _without_wildtype
                      sims_with_mutations=NA)


#A loop to go through all the sims from rep1 to rep32
for (i in 1:repfile) {
  # Load the file
  load(paste0("output/A_screen/SimResults_A_screen_rep_", i, ".RData"))
  # Create empty data frame for the values
  sims_summary <- data.frame(A = rep(A, each = l), #go through each file with l=1
                             PDtable = rep(i, length(A)),
                             entropy = NA,
                             entropy_without_wildtype=NA,
                             ed_average=NA,
                             ed_average_without_wildtype=NA,
                             bray_distance_average=NA,
                             bray_average_without_wildtype=NA,
                             sims_with_mutations=NA)
  # Analyze the loaded file using the analyseSims function with the empty data frame ready
  analysis_result <- analyseSims(sims)  
  
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "entropy"] <- analysis_result$entropy
  # Fill analysis_result$entropy_without_wildtype into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "entropy_without_wildtype"] <- analysis_result$entropy_without_wildtype
  # Fill analysis_result$ed_average into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "ed_average"] <- analysis_result$ed_average
  # Fill analysis_result$ed_average_without_wildtype into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "ed_average_without_wildtype"] <- analysis_result$ed_average_without_wildtype
  # Fill analysis_result$bray_average into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "bray_distance_average"] <- analysis_result$bray_distance_average
  # Fill analysis_result$bray_average_without_wildtype into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "bray_average_without_wildtype"] <- analysis_result$bray_average_without_wildtype
  # Fill analysis_result$sims_with_mutations into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "sims_with_mutations"] <- analysis_result$sims_with_mutations
}


#plot entropys
plot1<-ggplot(summary, aes(x = as.factor(A), y = entropy)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "A(Antibiotics concentration)", y = "Entropy") +
  theme_bw()
#plot ed_average
plot2<-ggplot(summary, aes(x = as.factor(A), y = ed_average)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "A(Antibiotics concentration)", y = "Eucliean distance") +
  theme_bw()
#plot entropy_without_wildtype
plot3<-ggplot(summary, aes(x = as.factor(A), y = entropy_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "A(Antibiotics concentration)", y = "Entropy") +
  theme_bw()
#plot ed_average_without_wildtype
plot4<-ggplot(summary, aes(x = as.factor(A), y = ed_average_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "A(Antibiotics concentration)", y = "Euclidean Distance") +
  theme_bw()

#plot bray_average
plot5<-ggplot(summary, aes(x = as.factor(A), y = bray_distance_average)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "A(Antibiotics concentration)", y = "Bray-Curtis Distance") +
  theme_bw()
#plot bray_average_without_wildtype
plot6<-ggplot(summary, aes(x = as.factor(A), y = bray_average_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "A(Antibiotics concentration)", y = "Bray-Curtis Distance ") +
  theme_bw()
#plot mutation sims_with_mutations  frequency 
plot7<-ggplot(summary, aes(x = as.factor(A), y = sims_with_mutations )) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "A(Antibiotics concentration)", y = "P(mutation spread)") +
  theme_bw()


#save plots
ggsave("plots/A_screen/entropy_A_screen.pdf", plot1, width =10, height = 10)
ggsave("plots/A_screen/ed_average_A_screen.pdf", plot2, width = 10, height = 10)
ggsave("plots/A_screen/entropy_without_wildtype_A_screen.pdf", plot3, width = 10, height = 10)
ggsave("plots/A_screen/ed_average_without_wildtype_A_screen.pdf", plot4, width = 10, height = 10)
ggsave("plots/A_screen/bray_average_A_screen.pdf", plot5, width = 10, height = 10)
ggsave("plots/A_screen/bray_average_without_wildtype_A_screen.pdf", plot6, width = 10, height = 10)
ggsave("plots/A_screen/sims_with_mutations_A_screen.pdf", plot7, width = 10, height = 10)

# save the database
#compount figure 4 A B C D without Wildtype 
save(summary, file = "output/A_screen/summary_A.RData")





library(cowplot)
#Plot all graphs contains mutants
plot_all_mutants<-plot_grid(plot7,plot3,plot4,plot6, labels = c("(a)", "(b)","(c)", "(d)"), hjust=0, vjust=1.3)
#Plot all graphs contains S and mutants
plot_all_S<-plot_grid(plot1,plot2,plot5, labels = c("(a)", "(b)","(c)"), hjust=0, vjust=1.3)

ggsave("plots/A_screen/plot_all_mutants.pdf", plot_all_mutants, width = 8, height = 8)
ggsave("plots/A_screen/plot_all_S.pdf", plot_all_S, width = 8, height=8)





