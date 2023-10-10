library(tidyverse)

#reset up all the useful parameters just in case it is other values
l <-1
n <- 2^(0:7)
m <- 1000
repfile<-32


# Load the necessary functions
source("scripts/NEW_ALL_functions.R")


# Create empty data frame for both entropy values and average_ed values
summary <- data.frame(n = rep(n, each = repfile),
                      PDtable = rep(1:repfile, length(n)),
                      entropy = NA,
                      entropy_without_wildtype=NA,
                      ed_average=NA,
                      ed_average_without_wildtype=NA,
                      bray_average=NA,
                      bray_average_without_wildtype=NA,
                      sims_with_mutations=NA)


#go through all the sims from rep1 to rep32
for (i in 1:repfile) {
  # Load the file
  load(paste0("SimResults_n_screen_rep_", i, ".RData"))
  # Create empty data frame
  sims_summary <- data.frame(n = rep(n, each = l),
                             PDtable = rep(i, length(n)),
                             entropy = NA,
                             entropy_without_wildtype=NA,
                             ed_average=NA,
                             ed_average_without_wildtype=NA,
                             bray_average=NA,
                             bray_average_without_wildtype=NA,
                             sims_with_mutations=NA)
  # Analyze the loaded file using the analyseSims function
  analysis_result <- analyseSims(sims)  
  
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "entropy"] <- analysis_result$entropy
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "entropy_without_wildtype"] <- analysis_result$entropy_without_wildtype
  # Fill analysis_result$ed_average into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "ed_average"] <- analysis_result$ed_average
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "ed_average_without_wildtype"] <- analysis_result$ed_average_without_wildtype
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "bray_average"] <- analysis_result$bray_average
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "bray_average_without_wildtype"] <- analysis_result$bray_average_without_wildtype
  # Fill analysis_result$entropy into the summary that has PDtable value equal to i
  summary[summary$PDtable == i, "sims_with_mutations"] <- analysis_result$sims_with_mutations
}


#plot entropys
plot1<-ggplot(summary, aes(x = as.factor(n), y = entropy)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "n(possible mutations)", y = "Entropy") +
  theme_bw()
#plot ed_average
plot2<-ggplot(summary, aes(x = as.factor(n), y = ed_average)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "n(possible mutations)", y = "Euclidean Distance") +
  theme_bw()
#plot entropy_without_wildtype
plot3<-ggplot(summary, aes(x = as.factor(n), y = entropy_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "n(possible mutations)", y = "Entropy ") +
  theme_bw()
#plot ed_average_without_wildtype
plot4<-ggplot(summary, aes(x = as.factor(n), y = ed_average_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "n(possible mutations)", y = "Euclidean Distance ") +
  theme_bw()
#plot bray_average_without_wildtype
plot5<-ggplot(summary, aes(x = as.factor(n), y = bray_average)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "n(possible mutations)", y = "Bray–Curtis Distance") +
  theme_bw()
#plot bray_average_without_wildtype
plot6<-ggplot(summary, aes(x = as.factor(n), y = bray_average_without_wildtype)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "n(possible mutations)", y = "Bray–Curtis Distance") +
  theme_bw()
#plot mutation sims_with_mutations  frequency 
plot7<-ggplot(summary, aes(x = as.factor(n), y = sims_with_mutations )) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2, height = 0, col = "blue", alpha = 0.5) +
  labs(x = "n(possible mutations)", y = "P(mutation spread)") +
  ylim(0, 1)+
  theme_bw()

#save plots
ggsave("plots/entropy_nscreen.pdf", plot1, width = 8, height = 8)
ggsave("plots/ed_average_nscreen.pdf", plot2, width = 8, height = 8)
ggsave("plots/entropy_without_wildtype_nscreen.pdf", plot3, width = 8, height = 8)
ggsave("plots/ed_average_without_wildtype_nscreen.pdf", plot4, width = 8, height = 8)
ggsave("plots/bray_average_nscreen.pdf", plot5, width = 8, height = 8)
ggsave("plots/bray_average_without_wildtype_nscreen.pdf", plot6, width = 8, height = 8)
ggsave("plots/sims_with_mutations_nscreen_10^6.pdf", plot7, width = 8, height = 8)

# save the database
save(summary, file = "output/summary_2.RData")

#combine plots of the mutants

library(cowplot)
#Plot all graphs contains mutants
plot_all_mutants<-plot_grid(plot7,plot3,plot4,plot6, labels = c("(a)", "(b)","(c)", "(d)"), hjust=0, vjust=1.3)
#Plot all graphs contains S and mutants
plot_all_S<-plot_grid(plot1,plot2,plot5, labels = c("(a)", "(b)","(c)"), hjust=0, vjust=1.3)

ggsave("plots/plot_all_mutants.pdf", plot_all_mutants, width = 8, height = 8)
ggsave("plots/plot_all_S.pdf", plot_all_S, width = 8, height=8)
