####----Before running ensure that in dataset the gene column is named as 'Gene', case sensitive.----------------------
# necessary library
library(tidyverse)

# Reading the dataset
data <- read.csv("filtered_genes_with_extra.csv")

# Defining function to calculate res_val(fi/bi)
calculate_res_val <- function(vec, m) {
    t <- sum(vec)
    b_i <- 1 / m
    f_i <- vec / t
    
    # Calculate res_val
    res_val <- f_i / b_i
    
    return(res_val)
}

# excluding gene column
m <- ncol(data) - 1

# Initializing vector to store column names
selected_high_columns <- c()
selected_low_columns <- c()

# Looping through each row (gene)
for (i in 1:nrow(data)) {
    gene_name <- make.names(data$Gene[i])  # Ensures valid file names for the gene
    values <- as.numeric(data[i, -1])  # Extracts the numeric expression values for the gene
    res_val <- calculate_res_val(values, m)  # Calculates res_val(fi/bi)
    
    sdv <- sd(res_val, na.rm = TRUE)  # Calculating Standard deviation, ignoring NA values
    mean_val <- mean(res_val, na.rm = TRUE)  # Calculating Mean, ignoring NA values
    
    # Finding values exceeding cutoff criteria to be considered overexpressed or under expressed i.e fi/bi > mean+sd or fi/bi < mean+sd
    high_values <- which(res_val > (mean_val + sdv))
    low_values <- which(res_val < (mean_val - sdv))
    
    #####-----------Separating out the genes as underexpressed or over expressed based on the nuber of profiles which are in high or low, there is cutoff ( how many profiles should be in high to be considered as overexpressed, change the cutoff as required and according to the dataset.------------########
    #### Storing column names if at least 40 values fulfill the condition (the cutoff here is 40, used for amr dataset2 where there are 196 profiles)
    if (length(high_values) > 40) {
        selected_high_columns <- c(selected_high_columns, gene_name)
    }

    if (length(low_values) > 40) {
        selected_low_columns <- c(selected_low_columns, gene_name)
    }
}

# Printing the overexpressed and underexpressed genes
print(selected_high_columns)
print(selected_low_columns)
#saving the list of genes which are overexpressed or underexpressed in separate files # can change the name of files as required
writeLines(paste(selected_high_columns, collapse = ","), "selected_high_genes.txt")
writeLines(paste(selected_low_columns, collapse = ","), "selected_low_genes.txt")
