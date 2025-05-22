# necessary libraries
library(tidyverse)
library(stats)

# Reading the data from amr dataset1 CSV
data <- read_csv("Japan_Lab_Transcriptomic_Data1.csv")

# Removing parent columns
#data <- data %>%
#  select(-c(Parent, `Parent (2nd)`))


# G-test (WRE) calculation function
calculate_g_test <- function(vec, m) {
  t <- sum(vec)  # Total sum of expression values
  b_i <- 1 / m  # Standard frequency in uniform distribution
  f_i <- vec / t  # Frequency of each expression value
  
  g_value <- 2 * sum(vec * log(f_i / b_i ))  # Adding a small value to avoid log(0)
  return(g_value)
}

# Finding total number of profiles (m)
m <- ncol(data) - 1  # Excluding the gene column

# Calculate G-test values for each gene in each sample
g_test_results <- data %>%
  rowwise() %>%
  mutate(
    G_value = calculate_g_test(c_across(-Gene), m)
  )

# Printing the G-test results
print(head(g_test_results))

# Saving the G-test results to a CSV file
write.csv(g_test_results, "G_test_results.csv", row.names = FALSE)
