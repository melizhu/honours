library(tidyverse)
i <- 2
j <- 888
singlesim <- sims[[i]]$results[[j]] |>
  as.tibble() |>
  pivot_longer(-time)
ggplot(singlesim) +
  geom_line(aes(x = time, y = value, colour = name))

