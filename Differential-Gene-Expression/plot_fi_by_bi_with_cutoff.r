#required libraries
library(tidyverse)

# Reading the dataset
data <- read.csv("filtered_genes_with_extra.csv")

# Remove parent columns if exist # comment out if not required
data <- data %>%
    select(-contains("Parent"))

#func to calculate res_val(fi/bi)
calculate_res_val <- function(vec, m) {
    t <- sum(vec)
    b_i <- 1 / m
    f_i <- vec / t
    
    # Calculating res_val
    res_val <- f_i / b_i
    
    return(res_val)
}

# Finding Number of samples (excluding gene column)
m <- ncol(data) - 1

# For loop to calculate res_val (fi/bi) for each gene and create histograms
for (i in 1:nrow(data)) {
    gene_name <- make.names(data$Gene[i])  # Making sure valid file names for the gene
    values <- as.numeric(data[i, -1])  # Extracting the numeric values for the gene (excluding the gene name)
    res_val <- calculate_res_val(values, m)  # Calculating res_val(fi/bi)
    sdv <- sd(res_val)  # Calculating standard deviation
    mean <- mean(res_val) # Calculating mean
    
    # Ploting the histogram
    plot <- ggplot(data.frame(res_val = res_val), aes(x = res_val-1)) +
        geom_histogram(bins = 30, fill = "skyblue", color = "black") +
        geom_vline(xintercept = mean - sdv, linetype = "dotted", color = "red", size = 1) +  # Adding left cutoff line
        geom_vline(xintercept = mean + sdv, linetype = "dotted", color = "red", size = 1) +  # Adding right cutoff line
        labs(
            title = paste("Histogram of res_val for", data$Gene[i]),
            x = "res_val",
            y = "Frequency"
        ) +
        theme_minimal()
    
    # Saving histogram to a file
    ggsave(
        filename = paste0(gene_name, "_histogram.png"),
        plot = plot,
        width = 8,
        height = 6
    )
}
