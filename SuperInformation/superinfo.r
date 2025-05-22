#necessary library
library(tidyverse)

#func to discretize the data
discretize_data <- function(seq, num_bins) {
    unique_vals <- length(unique(seq))
    if (unique_vals == 1) {
        return(rep(1, length(seq)))  # Return a single bin if all values are the same
    }
    breaks <- quantile(seq, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
    breaks <- unique(breaks)  # Ensure breaks are unique
    if (length(breaks) < num_bins + 1) {
        num_bins <- length(breaks) - 1  # Adjusting the number of bins
    }
    discretized_seq <- cut(seq, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    return(discretized_seq)
}

#func to calculate entropy
calculate_entropy <- function(seq) {
    # Computing probabilities
    probs <- table(seq) / length(seq)
    
    # calculating the entropy
    entropy <- -sum(probs * log2(probs))
    return(entropy)
}

# new discretize_data function
discretize_data <- function(seq, num_bins) {
    unique_vals <- length(unique(seq))
    if (unique_vals == 1) {
        return(rep(1, length(seq)))  # Return a single bin if all values are the same
    }
    breaks <- quantile(seq, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
    breaks <- unique(breaks)  # Ensure breaks are unique
    if (length(breaks) < num_bins + 1) {
        num_bins <- length(breaks) - 1  # Adjust number of bins
    }
    discretized_seq <- cut(seq, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    return(discretized_seq)
}

# Func to compute superinformation for each chunk
compute_chunk_superinformation <- function(sequence_chunk, block_size, num_bins) {
    num_blocks <- floor(length(sequence_chunk) / block_size)  # Number of complete blocks
    entropies <- numeric(num_blocks)  # Store entropy values
    
    # Computing entropy for each sub-block
    for (i in 1:num_blocks) {
        block <- sequence_chunk[((i - 1) * block_size + 1):(i * block_size)]  # Extract block
        entropies[i] <- calculate_entropy(block)  # Compute entropy
    }
    
    # Binning the entropy values
    entropy_bins <- discretize_data(entropies, num_bins)
    
    # Computing probabilities of entropy bins
    probs <- table(entropy_bins) / length(entropy_bins)
    
    # Computing entropy of entropy bins (superinformation)
    superinfo <- -sum(probs * log2(probs))
    
    return(superinfo)
}

# Main function to compute average superinformation across chunks
compute_superinformation <- function(sequence, chunk_size, block_size, num_bins) {
    num_chunks <- floor(length(sequence) / chunk_size)  # Number of complete chunks
    superinfo_values <- numeric(num_chunks)  # Store superinformation values for each chunk
    
    # Computing superinformation for each chunk
    for (i in 1:num_chunks) {
        chunk <- sequence[((i - 1) * chunk_size + 1):(i * chunk_size)]  # Extract chunk
        superinfo_values[i] <- compute_chunk_superinformation(chunk, block_size, num_bins)  # Compute superinformation
    }
    
    # Compute average superinformation across all chunks
    avg_superinfo <- mean(superinfo_values)
    
    return(avg_superinfo)
}

# can put another seq here
# Testing the functions with a sequence of 10000 'A's
sequence_as_characters <- rep("A", 10000)

#setting parameters
chunk_size <- 10000  # Setting chunk size
block_size <- 100  # Setting block size
num_bins <- 10  # Setting number of bins for discretization

# Computing average superinformation
avg_superinfo <- compute_superinformation(sequence_as_characters, chunk_size, block_size, num_bins)
print(avg_superinfo)
