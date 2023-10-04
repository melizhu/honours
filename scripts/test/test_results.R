
library(vegan)
library(readxl)


bray_curtis_dist <- function(x1, x2) {
  #Using the definition from Wikipedia.
  1 - (2 * sum(mapply(min, x1, x2))) / (sum(x1) + sum(x2))
}



custom_dist <- function(x, method, diag = FALSE, upper = FALSE) {
  as.dist(sapply(1:nrow(x), function(row1) {
    sapply(1:nrow(x), function(row2) {
      method(x[row1,], x[row2,])
    })
  }), diag = diag, upper = upper)
}

# Replace "sample_data.xlsx" with the actual path to your Excel file
my_data <- read_excel("test.xlsx")
# Assuming your features start from the second column, change the column index as needed.
distance <- vegdist(my_data[, -1], method = "euclidean") # Euclidean distance
bray_dissimilarity <- custom_dist(my_data[, -1], method = bray_curtis_dist)   # Bray-Curtis 
average_euclidean_distance <- mean(as.vector(distance), na.rm = TRUE)
average_bray_distance <- mean(as.vector(bray_distance), na.rm = TRUE)

