# necessary libraries
library(tidyverse)
library(ggplot2)
library(reshape2)

# Reading the data from CSV
data <- read.csv("filtered_genes.csv", header= TRUE, check.names = FALSE)  # Replace with your actual file name

# Function to calculate average of quadruplets
average_quadruplets <- function(df, column_prefix) {
    df %>%
        select(starts_with(column_prefix)) %>%
        rowMeans(na.rm = TRUE) %>%
        as.numeric()
}

# Calculate average for each sample and remove original columns
 data <- data %>%
     mutate(
         ATE_avg = average_quadruplets(., "3-ATE"),
         FOAE_avg = average_quadruplets(., "5-FOAE"),
         FUE_avg = average_quadruplets(., "5-FUE"),
         MPE_avg = average_quadruplets(., "6-MPE"),
         ABUE_avg = average_quadruplets(., "ABUE"),
         AFE_avg = average_quadruplets(., "AFE"),
         ATPE_avg = average_quadruplets(., "ATPE"),
         AZTE_avg = average_quadruplets(., "AZTE"),
         B_Cl_AlaE_avg = average_quadruplets(., "B-Cl-AlaE"),
         BSDE_avg = average_quadruplets(., "BSDE"),
         BZE_avg = average_quadruplets(., "BZE"),
         CBPCE_avg = average_quadruplets(., "CBPCE"),
         CCCPE_avg = average_quadruplets(., "CCCPE"),
         CMZE_avg = average_quadruplets(., "CMZE"),
         CPE_avg = average_quadruplets(., "CPE"),
         DCSE_avg = average_quadruplets(., "DCSE"),
         DVALE_avg = average_quadruplets(., "DVALE"),
         EDTAE_avg = average_quadruplets(., "EDTAE"),
         EME_avg = average_quadruplets(., "EME"),
         FOSE_avg = average_quadruplets(., "FOSE"),
         FTDE_avg = average_quadruplets(., "FTDE"),
         GAHE_avg = average_quadruplets(., "GAHE"),
         H2O2E_avg = average_quadruplets(., "H2O2E"),
         HSEE_avg = average_quadruplets(., "HSEE"),
         KME_avg = average_quadruplets(., "KME"),
         KTeE_avg = average_quadruplets(., "KTeE"),
         LVALE_avg = average_quadruplets(., "LVALE"),
         MECE_avg = average_quadruplets(., "MECE"),
         MMCE_avg = average_quadruplets(., "MMCE"),
         NFLXE_avg = average_quadruplets(., "NFLXE"),
         NiClE_avg = average_quadruplets(., "NiClE"),
         NITE_avg = average_quadruplets(., "NITE"),
         NMNOE_avg = average_quadruplets(., "NMNOE"),
         NQOE_avg = average_quadruplets(., "NQOE"),
         NVAE_avg = average_quadruplets(., "NVAE"),
         PHENE_avg = average_quadruplets(., "PHENE"),
         PLME_avg = average_quadruplets(., "PLME"),
         PMZE_avg = average_quadruplets(., "PMZE"),
         PSE_avg = average_quadruplets(., "PSE"),
         PURE_avg = average_quadruplets(., "PURE"),
         RFPE_avg = average_quadruplets(., "RFPE"),
         SDCE_avg = average_quadruplets(., "SDCE"),
         SHXE_avg = average_quadruplets(., "SHXE"),
         SSE_avg = average_quadruplets(., "SSE"),
         SXZE_avg = average_quadruplets(., "SXZE"),
         TETE_avg = average_quadruplets(., "TETE"),
         VCME_avg = average_quadruplets(., "VCME"),
         M9E_avg = average_quadruplets(., "M9E"),
         MDS42_avg = average_quadruplets(., "MDS42-n")
     )
 
 
 # Selecting only the relevant columns for MI
 selected_data <- data %>%
     select(ATE_avg, FOAE_avg, FUE_avg, MPE_avg, ABUE_avg, AFE_avg, ATPE_avg, AZTE_avg, BZE_avg, CBPCE_avg, CCCPE_avg, CMZE_avg, CPE_avg, DCSE_avg, DVALE_avg, EDTAE_avg, EME_avg, FOSE_avg, FTDE_avg, GAHE_avg, H2O2E_avg, HSEE_avg, KME_avg, KTeE_avg, LVALE_avg, MECE_avg, MMCE_avg, NFLXE_avg, NiClE_avg, NITE_avg, NMNOE_avg, NQOE_avg, NVAE_avg, PHENE_avg, PLME_avg, PMZE_avg, PSE_avg, PURE_avg, RFPE_avg, SDCE_avg, SHXE_avg, SSE_avg, SXZE_avg, TETE_avg, VCME_avg, M9E_avg, MDS42_avg)
  

# Function to calculate Mutual Information
mutual_information <- function(x, y) {
    # Ensuring only numeric data
    x <- as.numeric(x)
    y <- as.numeric(y)
    
    # Discretizing the data 
    num_bins <- 20  # Adjust the number of bins as needed
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

# Removing The gene name column
selected_data <- data %>% select(-gene)

# Calculating pairwise mutual information between all samples
mi_matrix <- matrix(0, nrow = ncol(selected_data), ncol = ncol(selected_data), 
                    dimnames = list(colnames(selected_data), colnames(selected_data)))

for (i in 1:ncol(selected_data)) {
    for (j in 1:ncol(selected_data)) {
        if (i == j) {
            mi_matrix[i, j] <- Inf  # Set mutual information between the same pair to infinity
        } else {
            mi_matrix[i, j] <- mutual_information(selected_data[[i]], selected_data[[j]])
        }
    }
}

# Melting the matrix for heatmap
melted_mi <- melt(mi_matrix)
colnames(melted_mi) <- c("Var1", "Var2", "MI")

# Saving the mutual information values to a CSV file, # change the name if required
write.csv(melted_mi, "filtered_mutual_information_values.csv", row.names = FALSE)


# Creating a heatmap with values
ggplot(melted_mi, aes(x = Var1, y = Var2, fill = MI)) + 
    geom_tile() + 
    geom_text(aes(label = ifelse(is.infinite(MI), "Inf", round(MI, 2))), size = 2, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = "Mutual Information Heatmap", x = "Sample", y = "Sample")

# Saving the heatmap to a file
ggsave("filtered_mutual_information_heatmap.png", width = 15, height = 10)
