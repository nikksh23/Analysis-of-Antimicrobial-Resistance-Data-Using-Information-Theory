#Required library
library(tidyverse)

# Read the dataset assumed that the file has one "Gene" column.
data <- read.csv("human_54samples_TMM_CPM.csv")

# Func to compute "fi/bi" for each gene samples,
calculate_res_val <- function(vec, m) {
t <- sum(vec)           # total expression for that gene
b_i <- 1 / m            # baseline value if expression were equally distributed
f_i <- vec / t          # fraction of total expression per sample
res_val <- f_i / b_i    # ratio: >1 means higher than baseline, <1 means lower
return(res_val)
}

# Number of samples (excluding the Gene column)
m <- ncol(data) - 1

# Creating a new output data frame to store for each gene:
# - The gene name,
# - The concatenated names of profiles (columns) that are overexpressed,
# - The concatenated names of profiles that are underexpressed.

result_df <- data.frame(
Gene = character(),
OverexpressedProfiles = character(),
UnderexpressedProfiles = character(),
stringsAsFactors = FALSE
)

# Iterating through each row (each gene) using for loop
for (i in 1:nrow(data)) {
# Using make.names in case gene names contain any spaces or special characters.
gene_name <- make.names(data$Gene[i])

# Extracting the numeric expression values for the gene, excluding the first column.
values <- as.numeric(data[i, -1])
# Calculate the fi/bi values for the profiles for this gene, If the expression were equally distributed, every profile would have res_val = 1.
res_val <- calculate_res_val(values, m)
# Compute the mean and standard deviation of the fi/bi values for this gene.
mean_val <- mean(res_val, na.rm = TRUE)
sd_val <- sd(res_val, na.rm = TRUE)

# Identifing which profiles exceed mean + sd (overexpressed)
high_indices <- which(res_val > (mean_val + sd_val))

# Identifing which profiles fall below mean - sd (underexpressed)
low_indices <- which(res_val < (mean_val - sd_val))

# Getting the sampl names, excluding the first "Gene" column.
profile_names <- colnames(data)[-1]

# Extractting the names of the profiles that meet the criterion.
high_profiles <- paste(profile_names[high_indices], collapse = ",")
low_profiles <- paste(profile_names[low_indices], collapse = ",")

# Defining the condition to classify gene e.g., "if at least X profiles are overexpressed, then record the gene". Appending the information for this gene to result dataframe.
result_df <- rbind(result_df, data.frame(
Gene = gene_name,
OverexpressedProfiles = high_profiles,
UnderexpressedProfiles = low_profiles,
stringsAsFactors = FALSE
))
}
# Exporting the results to a CSV file
write.csv(result_df, "gene_expression_profiles_comparison.csv", row.names = FALSE)