# necessary libraries
library(tidyverse)
library(ggplot2)
library(reshape2)

# Reading the data from CSV
data <- read_csv("filtered_file_newcutoff.csv")

# Function to calculate average of quadruplets
average_quadruplets <- function(df, column_prefix) {
  df %>%
    select(starts_with(column_prefix)) %>%
    rowMeans() %>%
    as.numeric()
}

# Calculating average for each sample
data <- data %>%
  mutate(
    CPZ_avg = average_quadruplets(., "CPZ"),
    CFIX_avg = average_quadruplets(., "CFIX"),
    AMK_avg = average_quadruplets(., "AMK"),
    NM_avg = average_quadruplets(., "NM"),
    DOXY_avg = average_quadruplets(., "DOXY"),
    CP_avg = average_quadruplets(., "CP"),
    AZM_avg = average_quadruplets(., "AZM"),
    TP_avg = average_quadruplets(., "TP"),
    ENX_avg = average_quadruplets(., "ENX"),
    CPFX_avg = average_quadruplets(., "CPFX")
  )

# Calculate average for parents # Only use when there's parent in the dataset otherwisw comment out
data <- data %>%
  mutate(
    Parent_avg = (Parent + `Parent (2nd)`) / 2
  )

# Select relevant columns for mutual information calculation # Remove Parent_avg if parent is not in dataset
selected_data <- data %>%
  select(Parent_avg, CPZ_avg, CFIX_avg, AMK_avg, NM_avg, DOXY_avg, CP_avg, AZM_avg, TP_avg, ENX_avg, CPFX_avg)

# Function to calculate Mutual Information
mutual_information <- function(x, y) {
  # Ensuring only numeric data
  x <- as.numeric(x)
  y <- as.numeric(y)

  # Discretizing the data
  num_bins <- 20  # Adjust the bin size here
  x_bins <- cut(x, breaks = num_bins, labels = FALSE)
  y_bins <- cut(y, breaks = num_bins, labels = FALSE)

  # Calculating joint and marginal probabilities
  joint_prob <- table(x_bins, y_bins) / length(x)
  prob_x <- prop.table(table(x_bins))
  prob_y <- prop.table(table(y_bins))

  # Adding small epsilon to avoid log2(0)
  epsilon <- .Machine$double.eps
  joint_prob <- joint_prob + epsilon
  prob_x <- prob_x + epsilon
  prob_y <- prob_y + epsilon

  # Calculating mutual information
  mi <- sum(joint_prob * log2(joint_prob / (outer(prob_x, prob_y))))

  return(mi)
}

# Calculating pairwise mutual information between all samples
mi_matrix <- matrix(0, nrow = ncol(selected_data), ncol = ncol(selected_data), 
                    dimnames = list(colnames(selected_data), colnames(selected_data)))

for (i in 1:ncol(selected_data)) {
  for (j in 1:ncol(selected_data)) {
    mi_matrix[i, j] <- mutual_information(selected_data[[i]], selected_data[[j]])
  }
}

# Melting the matrix for heatmap
melted_mi <- melt(mi_matrix)
colnames(melted_mi) <- c("Var1", "Var2", "MI")

# Saving the mutual information values to a CSV file # change the name here as desired
write.csv(melted_mi, "mutual_information_values_newcutoff.csv", row.names = FALSE)


# Creating a heatmap using ggplot, can change parameters as required
ggplot(melted_mi, aes(x = Var1, y = Var2, fill = MI)) + 
  geom_tile() + 
  geom_text(aes(label = round(MI, 2)), size = 2, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Mutual Information Heatmap", x = "Sample", y = "Sample")

# Save the heatmap to a file # can change the name of the heatmap saved as required
ggsave("mutual_information_heatmap_newcutoff.png", width = 15, height = 10)

