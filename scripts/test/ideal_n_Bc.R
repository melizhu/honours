#Calculate the theoretical maximum Bray_Curtis distance for different n

# Define the number of dimensions and points
num_dimensions <- 2^(0:7)
num_points <- num_dimensions

for (k in 1:length(num_points)){
# Initialize an empty matrix to store points
points_matrix <- diag(1, nrow = num_points[k], ncol = num_dimensions[k])

# Initialize an empty matrix to store distances
distance_matrix <- matrix(NA, nrow = num_points[k], ncol = num_points[k])

# Calculate the Bray-Curtis distance between all combinations
for (i in 1:num_points[k]) {
  for (j in 1:num_points[k]) {
      numerator <- sum(abs(points_matrix[i,] - points_matrix[j,]))
      denominator <- sum(abs(points_matrix[i,] + points_matrix[j,]))
      distance_matrix[i, j] <- numerator / denominator
    
   }
 }
# Calculate the average distance
average_distance[k] <- sum(distance_matrix) / (num_points[k] * num_points[k] )
}
# Print the average distance
print(average_distance)
