## This is for amr dataset2####
#necessary library
library(tidyverse)
library(ggplot2)

# function to discretize data
discretize_data <- function(seq, num_bins) {
  breaks <- quantile(seq, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
  discretized_seq <- cut(seq, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  return(discretized_seq)
}

# Func to calculate entropy
calculate_entropy <- function(seq) {
  probs <- table(seq) / length(seq)
  entropy <- -sum(probs * log2(probs + .Machine$double.eps), na.rm = TRUE)
  return(entropy)
}

# Func to calculate joint entropy of two
calculate_joint_entropy <- function(seq1, seq2) {
  joint_table <- table(seq1, seq2) / length(seq1) + .Machine$double.eps
  joint_entropy <- -sum(joint_table * log2(joint_table), na.rm = TRUE)
  return(joint_entropy)
}

# Func to calculate joint entropy of three variables
calculate_triple_entropy <- function(seq1, seq2, seq3) {
  joint_table <- table(seq1, seq2, seq3) / length(seq1) + .Machine$double.eps
  joint_entropy <- -sum(joint_table * log2(joint_table), na.rm = TRUE)
  return(joint_entropy)
}

##### Func to calculate three-variable Mutual Information (MMI)######
calculate_mmi <- function(data, A, B, C) {
  # Extracting the expresseion values for A, B, and C ( A, B & C are any three profiles picked from the dataset)
  seqA <- data[[A]]
  seqB <- data[[B]]
  seqC <- data[[C]]
  
  # calling the func for individual entropies
  HA <- calculate_entropy(seqA)
  HB <- calculate_entropy(seqB)
  HC <- calculate_entropy(seqC)
  
  # Calling the func for pairwise joint entropies
  HAB <- calculate_joint_entropy(seqA, seqB)
  HBC <- calculate_joint_entropy(seqB, seqC)
  HCA <- calculate_joint_entropy(seqC, seqA)
  
  # Calling the func for joint entropy for all three variables
  HABC <- calculate_triple_entropy(seqA, seqB, seqC)
  
  # Calculating Pairwise MI and MMI
  MI12 <- HA + HB - HAB  # MI between seqA and seqB
  MI23 <- HB + HC - HBC  # MI between seqB and seqC
  MI13 <- HA + HC - HCA  # MI between seqA and seqC
  MMI <- HA + HB + HC - (HAB + HBC + HCA) + HABC
  
  return(list(MI12 = MI12, MI23 = MI23, MI13 = MI13, MMI = MMI))
}

data <- read.csv("filtered_genes.csv", header= TRUE, check.names = FALSE)

### Function to average of quadruplets####
average_quadruplets <- function(df, column_prefix) {
  df %>%
    select(starts_with(column_prefix)) %>%
    rowMeans() %>%
    as.numeric()
}

# Calculating average for each sample
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
 
 
 #Selecting averaged columns for MI
 selected_data <- data %>%
     select(ATE_avg, FOAE_avg, FUE_avg, MPE_avg, ABUE_avg, AFE_avg, ATPE_avg, AZTE_avg, BZE_avg, CBPCE_avg, CCCPE_avg, CMZE_avg, CPE_avg, DCSE_avg, DVALE_avg, EDTAE_avg, EME_avg, FOSE_avg, FTDE_avg, GAHE_avg, H2O2E_avg, HSEE_avg, KME_avg, KTeE_avg, LVALE_avg, MECE_avg, MMCE_avg, NFLXE_avg, NiClE_avg, NITE_avg, NMNOE_avg, NQOE_avg, NVAE_avg, PHENE_avg, PLME_avg, PMZE_avg, PSE_avg, PURE_avg, RFPE_avg, SDCE_avg, SHXE_avg, SSE_avg, SXZE_avg, TETE_avg, VCME_avg, M9E_avg, MDS42_avg)
  


# Calling function to Discretize the data
num_bins <- 20  # Set number of bins for discretization # can chamge according to need
discretized_data <- selected_data %>%
  mutate(across(everything(), ~ discretize_data(.x, num_bins)))

# Calculating MMI and Pairwise MI for unique combinations of three columns by calling the MMI function and iterating through the dataset
unique_combinations <- combn(colnames(discretized_data), 3, simplify = FALSE)  # Generate unique combinations
results <- tibble(A = character(), B = character(), C = character(), MI12 = numeric(), MI23 = numeric(), MI13 = numeric(), MMI = numeric())

#for loop for iterating and savin the results in results dataframe
for (combo in unique_combinations) {
  mi_values <- calculate_mmi(discretized_data, combo[1], combo[2], combo[3])
  results <- results %>% add_row(A = combo[1], B = combo[2], C = combo[3], MI12 = mi_values$MI12, MI23 = mi_values$MI23, MI13 = mi_values$MI13, MMI = mi_values$MMI)
}


#Saving results to a CSV file # can change the name as required
write.csv(results, "three_variable_mi_filtered_genes.csv", row.names = FALSE)


###--- The following setps are for plotting which is not necessarily required and can be skipped ------####
# Figuring out combination identifier for plotting
results <- results %>%
    mutate(Combination = paste(A, B, C, sep = "_"))

# using Bar plot for MMI values
ggplot(results, aes(x = Combination, y = MMI)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
    labs(title = "Three-Variable Mutual Information (MMI)", 
         x = "Variable Combinations", 
         y = "MMI Value")


