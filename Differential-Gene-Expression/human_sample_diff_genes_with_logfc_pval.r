# necessary libraries
library(tidyverse)
library(stringr)

# Reading the dataset
data <- read.csv("human_54samples_TMM_CPM.csv")

# Func to compute res_val(fi/bi)
calculate_res_val <- function(vec, m) {
    t <- sum(vec)
    b_i <- 1 / m
    f_i <- vec / t
    res_val <- f_i / b_i
    return(res_val)
}

# Extracting sample names
sample_names <- colnames(data)[-1]

# Identifing unique sample names(not considering replicates) (removing ".1", ".2" replicate labels)
base_profiles <- unique(str_replace(sample_names, "\\.[0-9]+$", ""))

# Adding a small pseudo-count to avoid log2(0) errors
pseudo_count <- 1e-6

# Initializing list for results
results_list <- list()

# Iterating and identifying the differentially expressed genes
for(profile in base_profiles) {
  
  # Identifing replicates in the test group
  case_indices <- which(str_detect(sample_names, paste0("^", profile, "(\\.\\d+)?$")))
  # setting the rest of samples as control
  control_indices <- setdiff(seq_along(sample_names), case_indices)
  
  # Initializing lists
  over_expressed_genes <- data.frame(Gene = character(), Log2FC = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
  under_expressed_genes <- data.frame(Gene = character(), Log2FC = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

#Iterating through each of gene, using their expresion values to check if the text is over expressed or underexpressed compared to control
  for(i in 1:nrow(data)) {
    gene_name <- data$Gene[i]
    values <- as.numeric(data[i, -1])
    
    # Computing res_val(fi/bi) for all samples
    res_val <- calculate_res_val(values, length(sample_names))
    
    # Extracting control values and compute mean adn standart deviation for cutoffs
    control_val <- res_val[control_indices]
    control_mean <- mean(control_val, na.rm = TRUE)
    control_sd <- sd(control_val, na.rm = TRUE)
    
    # Computing average res_val for the test group replicates
    test_val <- res_val[case_indices]
    test_avg_res_val <- mean(test_val, na.rm = TRUE)
    
    # Applying pseudo-count if test_avg_res_val is zero to avoid log2(0) issues
    test_avg_res_val <- ifelse(test_avg_res_val == 0, pseudo_count, test_avg_res_val)
    control_mean <- ifelse(control_mean == 0, pseudo_count, control_mean)

    # Calculating log2 fold-change (LFC) 
    log2_fc <- if (!is.na(test_avg_res_val) && !is.na(control_mean)) {
      log2(test_avg_res_val / control_mean)
    } else {
      NA  # Marking as missing if calculations are invalid
    }
    
    # Performing Welch’s t-test between test and control groups
    if (length(test_val) > 1 && length(control_val) > 1) {
      t_test <- t.test(test_val, control_val, var.equal = FALSE)  # Welch's t-test using t-test function
      p_value <- t_test$p.value
    } else {
      p_value <- NA  # giving NA if insufficient data points
    }
    
    # Separating genes into overexpressed or underexpressed using there fi/bi values using cutoff mean+sd and mean-sd
    if(test_avg_res_val > (control_mean + control_sd)) {
      over_expressed_genes <- rbind(over_expressed_genes, data.frame(Gene = gene_name, Log2FC = log2_fc, P_Value = p_value))
    }
    if(test_avg_res_val < (control_mean - control_sd)) {
      under_expressed_genes <- rbind(under_expressed_genes, data.frame(Gene = gene_name, Log2FC = log2_fc, P_Value = p_value))
    }
  }
  
  # Storing results
  results_list[[profile]] <- list(
    profile = profile,
    over_expressed = over_expressed_genes,
    under_expressed = under_expressed_genes
  )
}

# Saving the results to CSV files and making sure there are genes, if not it will abort the saving process and say no genes found, overexpressed and underexpressed genes are stored in two different csv files along with their logFC and P values.
for(profile in names(results_list)) {
  if(nrow(results_list[[profile]]$over_expressed) > 0) {
    write.csv(results_list[[profile]]$over_expressed, paste0("overexpressed_genes_with_log2FC_WelchT_", profile, ".csv"), row.names = FALSE)
  } else {
    cat("No overexpressed genes found for", profile, "\n")
  }

  if(nrow(results_list[[profile]]$under_expressed) > 0) {
    write.csv(results_list[[profile]]$under_expressed, paste0("underexpressed_genes_with_log2FC_WelchT_", profile, ".csv"), row.names = FALSE)
  } else {
    cat("No underexpressed genes found for", profile, "\n")
  }
}

# Printing Summary of the results found
for(profile in names(results_list)) {
  cat("Profile:", profile, "\n")
  cat("  Overexpressed genes (n =", nrow(results_list[[profile]]$over_expressed), ")\n")
  cat("  Underexpressed genes (n =", nrow(results_list[[profile]]$under_expressed), ")\n\n")
}
